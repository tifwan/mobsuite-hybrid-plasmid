
# open sequence QC metadata
import pandas as pd

#  Find the species for each isolate from the QC metadata
metadata = pd.read_csv("/scratch/esnitkin_root/esnitkin1/tifwan/Project_MDHHS_genomics/2026-09-24_hybrid_mobsuite/mobsuite-hybrid-plasmid/meta/master_qc_summary.csv")
taxonomy = metadata[["Sample", "Species"]]

# rename columns to match expected format
taxonomy = taxonomy.rename(columns={"Sample": "id", "Species": "organism"})

# write taxonomy to tsv
taxonomy.to_csv("/scratch/esnitkin_root/esnitkin1/tifwan/Project_MDHHS_genomics/2026-09-24_hybrid_mobsuite/mobsuite-hybrid-plasmid/meta/taxonomy.tsv", sep="\t", index=False )