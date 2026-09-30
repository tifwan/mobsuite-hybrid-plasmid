# Author: Tiffany, Auden Bahr
# Date: Sep. 29, 2026

import os
import glob
from datetime import datetime
RUN_TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")



configfile: "config/config_cluster.yaml"
container: "docker://" + config["mob_image"]

OUTPUT_DIR = config["output_dir"]

# Discover isolates from what the upstream mob_recon/rename pipeline has
# already finished -- this script only needs to know which ones are ready
# for QC/rename/clustering, not re-derive the full assembly list.
ISOLATES = sorted(
    os.path.basename(os.path.dirname(p))
    for p in glob.glob(os.path.join(OUTPUT_DIR, "*", "*_renamed.done"))
)
isolate_prefix = "MI_KPC"


# ==============================================================================
rule all:
    input:
        OUTPUT_DIR + "/mob_cluster/database_check.txt", 


# ---- 1. QC FIRST, on the raw untouched mob_recon output --------------------
# Contig count only means anything while each plasmid_<ID>.fasta is still its
# own separate file -- check it here, before anything gets renamed or
# combined, so no group-tracking machinery is needed at all. Qualifying
# (single-contig) files are copied, untouched, into qc_passed/.
rule assembly_qc:
    input:
        renamed = OUTPUT_DIR + "/{isolate}/{isolate}_renamed.done",
    output:
        qc_dir = directory(OUTPUT_DIR + "/{isolate}/mob_cluster_prep/qc_passed"),
    params:
        recon_dir = OUTPUT_DIR + "/{isolate}/mob_recon",
    log:
        config["log_dir"] + "/{isolate}/assembly_qc_" + RUN_TIMESTAMP + ".log",
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p {output.qc_dir}
        shopt -s nullglob
        for f in {params.recon_dir}/plasmid_*.fasta; do
            n=$(grep -c "^>" "$f")
            if [ "$n" -eq 1 ]; then
                cp "$f" {output.qc_dir}/
            else
                echo "excluding $f: $n contigs (fragmented, not closed)"
            fi
        done
        shopt -u nullglob
        }} > {log} 2>&1
        """
 
 
# ---- 2. Rename headers on ONLY the QC-passed files, then combine -----------
# Since QC already ran on the untouched files, this step can freely combine
# everything for the isolate into one output without losing any information
# the QC step needed -- that information was already consumed in step 1.
#
# This is also where the host-taxonomy file gets built, one isolate at a
# time, because {params.isolate} is a known, exact string right here -- doing
# it later (after headers from many isolates are combined) would mean
# prefix-matching isolate names back out of strings like "MI_KPC_1_AB978_5",
# which is genuinely awkward since isolate names themselves contain
# underscores. Building it here instead is just one extra lookup + echo per
# record, no separate script needed.
#
# Assumptions to verify against your real files:
#  - mobtyper_results.txt's first column is the sequence ID, matching the
#    fasta header exactly (if not column 1, or headers carry more than just
#    an ID after '>', adjust the awk/sed below).
#  - taxonomy.tsv has a header row with columns named `id` and `organism`
#    (tab-separated), one row per isolate.
#  - mob_cluster's expected -t format is two tab-separated columns (sequence
#    id, organism) -- this hasn't been confirmed against `mob_cluster --help`
#    or an example file; adjust if it wants something else.
rule rename_fasta_headers:
    input:
        qc_dir       = OUTPUT_DIR + "/{isolate}/mob_cluster_prep/qc_passed",
        taxonomy_tsv = config["taxonomy_file"],
    output:
        fasta    = OUTPUT_DIR + "/{isolate}/mob_cluster_prep/{isolate}_qc_passed.fasta",
        report   = OUTPUT_DIR + "/{isolate}/mob_cluster_prep/{isolate}_qc_passed_report.txt",
        taxonomy = OUTPUT_DIR + "/{isolate}/mob_cluster_prep/{isolate}_taxonomy.tsv",
    params:
        recon_dir = OUTPUT_DIR + "/{isolate}/mob_recon",  # mobtyper_results.txt lives here
        isolate   = "{isolate}",
    log:
        config["log_dir"] + "/{isolate}/rename_headers_" + RUN_TIMESTAMP + ".log",
    shell:
        r"""
        {{
        set -euo pipefail
        : > {output.fasta}
        : > {output.taxonomy}
        head -n 1 {params.recon_dir}/mobtyper_results.txt > {output.report} 2>/dev/null || : > {output.report}
 
        # Look this isolate's organism up ONCE -- it's the same for every
        # plasmid belonging to it.
        organism=$(awk -F'\t' -v id="{params.isolate}" 'NR>1 && $1==id {{print $2; exit}}' {input.taxonomy_tsv})
        if [ -z "$organism" ]; then
            echo "WARNING: no taxonomy.tsv row found for isolate {params.isolate}" >&2
        fi
 
        shopt -s nullglob
for f in {input.qc_dir}/plasmid_*.fasta; do
    id=$(basename "$f" .fasta | sed 's/^plasmid_//')
    orig_header=$(grep "^>" "$f" | head -1 | sed 's/^>//')
    mobtyper_id="{params.isolate}:${{id}}"
    new_header="{params.isolate}_${{id}}_${{orig_header}}"

    sed "s/^>.*/>${{new_header}}/" "$f" >> {output.fasta}

    awk -F'\t' -v OFS='\t' -v old="$mobtyper_id" -v new="$new_header" \
        'NR>1 && $1==old {{ $1=new; print }}' {params.recon_dir}/mobtyper_results.txt >> {output.report}

    printf '%s\t%s\n' "$new_header" "$organism" >> {output.taxonomy}
done
shopt -u nullglob
        }} > {log} 2>&1
        """
 
 
# ---- 3. Aggregate every isolate's QC-passed plasmids into one batch --------
rule concatenate_qc_passed_plasmids:
    input:
        fastas     = expand(OUTPUT_DIR + "/{isolate}/mob_cluster_prep/{isolate}_qc_passed.fasta", isolate=ISOLATES),
        reports    = expand(OUTPUT_DIR + "/{isolate}/mob_cluster_prep/{isolate}_qc_passed_report.txt", isolate=ISOLATES),
        taxonomies = expand(OUTPUT_DIR + "/{isolate}/mob_cluster_prep/{isolate}_taxonomy.tsv", isolate=ISOLATES),
    output:
        fasta    = OUTPUT_DIR + "/mob_cluster_prep/all_new_plasmids.fasta",
        report   = OUTPUT_DIR + "/mob_cluster_prep/all_new_plasmids_mobtyper_report.txt",
        taxonomy = OUTPUT_DIR + "/mob_cluster_prep/all_new_plasmids_host_taxonomy.tsv",
    log:
        config["log_dir"] + "/concatenate_qc_passed_" + RUN_TIMESTAMP + ".log",
    shell:
        r"""
        {{
        set -euo pipefail
        : > {output.fasta}
        for f in {input.fastas}; do
            cat "$f" >> {output.fasta}
        done
 
        : > {output.report}
        header_written=0
        for r in {input.reports}; do
            if [ -s "$r" ]; then
                if [ "$header_written" -eq 0 ]; then
                    head -n 1 "$r" > {output.report}
                    header_written=1
                fi
                tail -n +2 "$r" >> {output.report}
            fi
        done
 
        : > {output.taxonomy}
        printf 'id\torganism\n' > {output.taxonomy}
        for t in {input.taxonomies}; do
        cat "$t" >> {output.taxonomy}
        done
        }} > {log} 2>&1
        """
 
 
 
# Update MOBsuite database by modifying files loaded by conda.
#rule mob_cluster:
#    input:
#        all_plasmids = all_plasmids_fasta,
#        new_plasmids = OUTPUT_DIR + "/mob_cluster_prep/all_new_plasmids.fasta",
#        mobtyper =  OUTPUT_DIR + "/mob_cluster_prep/all_new_plasmids_mobtyper_report.txt",
#        tfile = OUTPUT_DIR + "/mob_cluster_prep/all_new_plasmids_host_taxonomy.tsv",
#        old_clusters = str(config["mob_db"]) + "/clusters.txt",

#    output:
#        updated_db = OUTPUT_DIR + "/mob_cluster/database_check.txt",
#    params:
#        outdir = OUTPUT_DIR + "/mob_cluster",
#        old_db = str(config["mob_db"]),
#        copyfiles = str(config["mob_db"]) + "/*"

#    log:
#        "logs/mob_cluster/mob_cluster_" + RUN_TIMESTAMP + ".log"
#    shell:
#        """
#        mob_cluster --mode update -f {input.new_plasmids} -p {input.mobtyper} -t {input.tfile} -c {input.old_clusters} -r {input.all_plasmids} --outdir {params.outdir} >> {log} 2>&1
#        cp {params.outdir}/clusters.txt {params.old_db}/clusters.txt >> {log} 2>&1
#        cp {params.outdir}/references_updated.fasta {params.old_db}/ncbi_plasmid_full_seqs.fas >> {log} 2>&1
#        makeblastdb -in {params.old_db}/ncbi_plasmid_full_seqs.fas -dbtype nucl >> {log} 2>&1
#        mash sketch -i {params.old_db}/ncbi_plasmid_full_seqs.fas >> {log} 2>&1
#        grep "{isolate_prefix}" {params.old_db}/clusters.txt > {output.updated_db}
#        """

rule concat_all_plasmids_fasta:
    input: 
        old_plasmids_fasta = config["mob_db"] + "/ncbi_plasmid_full_seqs.fas",
        new_plasmids = OUTPUT_DIR + "/mob_cluster_prep/all_new_plasmids.fasta",
    output:
        concatenated = OUTPUT_DIR + "/mob_cluster_prep/all_plasmids.fasta",
    shell:
        """
        cat {input.old_plasmids_fasta} {input.new_plasmids} > {output.concatenated}
        """
        
# Update MOBsuite database by modifying files loaded by conda.
rule mob_cluster:
    input:
        all_plasmids = OUTPUT_DIR + "/mob_cluster_prep/all_plasmids.fasta",
        new_plasmids = OUTPUT_DIR + "/mob_cluster_prep/all_new_plasmids.fasta",
        mobtyper =  OUTPUT_DIR + "/mob_cluster_prep/all_new_plasmids_mobtyper_report.txt",
        tfile = OUTPUT_DIR + "/mob_cluster_prep/all_new_plasmids_host_taxonomy.tsv",
        old_clusters = str(config["mob_db"]) + "/clusters.txt",
    output:
        updated_db = OUTPUT_DIR + "/mob_cluster/database_check.txt",
    params:
        outdir = OUTPUT_DIR + "/mob_cluster",
        old_db = str(config["mob_db"]),
        copyfiles = str(config["mob_db"]) + "/*"

    log:
        "logs/mob_cluster/mob_cluster_" + RUN_TIMESTAMP + ".log"
    shell:
        """
        mob_cluster --mode update -f {input.new_plasmids} -p {input.mobtyper} -t {input.tfile} -c {input.old_clusters} -r {input.all_plasmids} --outdir {params.outdir} >> {log} 2>&1
        shopt -s nullglob
        for f in {params.outdir}/references_updated.fasta*; do
            new="${{f/references_updated.fasta/ncbi_plasmid_full_seqs.fas}}"
            mv "$f" "$new"
        done
        shopt -u nullglob
        makeblastdb -in {params.outdir}/ncbi_plasmid_full_seqs.fas -dbtype nucl >> {log} 2>&1
        cp {params.outdir}/ncbi* {params.outdir}/clusters.txt {params.old_db} >> {log} 2>&1
        grep "{isolate_prefix}" {params.old_db}/clusters.txt > {output.updated_db}
        """