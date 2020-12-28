
import pandas as pd
import glob


## --------------------------------------------------------------------------------
## global parameters from config file

PREFIX=config["prefix"] ## output prefix


## --------------------------------------------------------------------------------
## helpers

unit_df = pd.read_table(config["units"], comment="#").set_index(["sampleId"], drop=False)
SAMPLES = unit_df.index.unique()


## --------------------------------------------------------------------------------
## targets

rule all:
    input:
        "plots/" + PREFIX + ".hits.pdf",
        "tables/" + PREFIX + ".summary.tsv.gz"

          
## --------------------------------------------------------------------------------
## rules

rule summary_all:
    input:
        tsv=expand("tables/{sample}/" + PREFIX + ".summary.tsv.gz", sample = SAMPLES),
        stamp=expand("stages/{sample}.summary.done", sample = SAMPLES),
    output:
        pdf="plots/" + PREFIX + ".hits.pdf",
        tsv="tables/" + PREFIX + ".summary.tsv.gz"
    shell:
        """
        Rscript src/summaryAll.R {PREFIX} {input.tsv}
        """

