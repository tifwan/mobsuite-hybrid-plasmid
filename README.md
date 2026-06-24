# MOB-suite hybrid pipeline

Runs MOB-suite on **hybrid assemblies only**:

1. `mob_recon` per isolate (plasmid reconstruction + typing — `mobtyper_results.txt`
   is produced automatically).
2. Renames every `plasmid_*.fasta` and `chromosome.fasta` by prefixing the isolate
   name, e.g. `plasmid_AA337.fasta` → `IMPALA_1000_plasmid_AA337.fasta`.
3. Concatenates all typing reports into a single `output/mob_typer.txt`.

No Illumina, no reference-database update, no `mob_cluster`.

## Layout
```
config/
  config.yaml            # paths + settings (edit this)
  cluster.json           # SLURM resource config
  isolates_hybrid.csv    # one column "assembly_dir", one isolate per row
workflow/
  Snakefile
  envs/mob.yaml          # conda fallback
  scripts/
logs/                    # per-rule logs, timestamped (created on run)
output/
  mob_typer.txt
  <isolate>/
    mob_recon/...                       # raw mob_recon output
    <isolate>_chromosome.fasta          # renamed
    <isolate>_plasmid_*.fasta           # renamed
```

## Setup
1. Edit `config/config.yaml`:
   - `assembly_dir_hybrid` — where your hybrid assemblies live
   - `assembly_pattern` — how to find each FASTA (default `{isolate}.fasta`)
   - `chrom_filt_set` — optional chromosome filter; leave `""` to use the default
2. List your isolates in `config/isolates_hybrid.csv` under the `assembly_dir` column.

## Run

### With the Docker image via Singularity/Apptainer (recommended)
Snakemake pulls `kbessonov/mob_suite` and runs it through Singularity. Docker
itself needs root, so on HPC use Singularity/Apptainer:

```bash
snakemake -s workflow/Snakefile --use-singularity --cores 8
```

The MOB-suite databases (~large) download automatically on first run inside the
container. To persist them across runs, bind a host directory and point
MOB-suite at it (add `--singularity-args "-B /path/to/db:/db"` and pass
`-d /db` via the rule, or run `mob_init` once into a bound location).

### Conda fallback (no Singularity)
```bash
snakemake -s workflow/Snakefile --use-conda --cores 8
```
(Add `conda: "envs/mob.yaml"` to each rule, or it uses the container directive.)

### On a SLURM cluster
```bash
snakemake -s workflow/Snakefile --use-singularity \
  --cluster "sbatch -A {cluster.account} -p {cluster.partition} -t {cluster.time} \
             -N {cluster.nodes} -n {cluster.ntasks} -c {cluster.cpus-per-task} \
             --mem {cluster.mem} -J {cluster.job-name} \
             -o {cluster.output} -e {cluster.error}" \
  --cluster-config config/cluster.json \
  --jobs 20
```

### Dry run (check the plan first)
```bash
snakemake -s workflow/Snakefile -n
```

## Logs
Every rule writes a timestamped log under `logs/<isolate>/`. Cluster scheduler
stdout/stderr go to `logs/cluster/`. Check these first when troubleshooting.

## Notes
- `assembly_pattern` lets you match either a flat layout (`{isolate}.fasta`) or
  the original nested one (`{isolate}/{isolate}_flye_medaka_polypolish.fasta`).
- All plasmids are renamed (no single-contig/circular filtering). If you only
  want complete plasmids, add that filter back into `rename_outputs`.
- Swap the MOB-suite version via `mob_image` in the config (e.g. `:3.0.3`).
