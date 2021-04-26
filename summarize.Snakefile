
import pandas as pd
import glob


## --------------------------------------------------------------------------------
## global parameters from config file

PREFIX=config["prefix"] ## output prefix
TARGETS_PRIORITY=config["targets_priority"] ## list of priority genera


## --------------------------------------------------------------------------------
## helpers

unit_df = pd.read_table(config["units"], comment="#").set_index(["sampleId"], drop=False)
SAMPLES = unit_df.index.unique()


## --------------------------------------------------------------------------------
## targets

rule all:
    input:
        "plots/" + PREFIX + ".targets_priority.matrix.pdf",
        "plots/" + PREFIX + ".targets_priority.pdf",
        "tables/" + PREFIX + ".summary.tsv.gz",
        "tables/" + PREFIX + ".targets_priority.summary.tsv.gz"

          
## --------------------------------------------------------------------------------
## rules

rule summary_all:
    input:
        tsv=expand("tables/{sample}/" + PREFIX + ".summary.tsv.gz", sample = SAMPLES),
        stamp=expand("stages/{sample}.summary.done", sample = SAMPLES),
    output:
        "plots/" + PREFIX + ".targets_priority.matrix.pdf",
        "plots/" + PREFIX + ".targets_priority.pdf",
        "tables/" + PREFIX + ".summary.tsv.gz",
        "tables/" + PREFIX + ".targets_priority.summary.tsv.gz"
    shell:
        """
        Rscript src/summaryAll.R {PREFIX} {TARGETS_PRIORITY} {input.tsv}
        """

