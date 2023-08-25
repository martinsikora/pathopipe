import pandas as pd
import glob


## --------------------------------------------------------------------------------
## global parameters from config file

PREFIX = config["prefix"]  ## output prefix
DB = config["db"]  ## database root dir
GENFILE = config["genera"]  ## target genera file
KMERS = config["kmers"]  ## minimum kmers for target assemblies
MQ = config["MQ"]  ## MQ cutoff for mapped BAMs
TMP_DIR = config["tmp_dir"]  ## path to location for temporary files
MM2_PARAM = config["mm2_param"]
METADMG_CPP = config["metaDMG-cpp"]  ## path to metaDMG cpp module


## --------------------------------------------------------------------------------
## helpers

unit_df = pd.read_table(config["units"], comment="#").set_index(
    ["sampleId"], drop=False
)
SAMPLES = unit_df.index.unique()

db_df = pd.read_table(DB + "/library.seqInfo.tsv").set_index(["assemblyId"], drop=False)

gen_df = pd.read_table(GENFILE, comment="#", header=None)
GENERA = gen_df[[1]]


## --------------------------------------------------------------------------------
## functions


def get_mmi(wildcards):
    return DB + "/" + db_df.loc[(wildcards.assembly), "mmi"]


def get_fa(wildcards):
    return DB + "/" + db_df.loc[(wildcards.assembly), "fa"]


def get_fq(wildcards):
    return unit_df.loc[(wildcards.sample), "fq"]


## --------------------------------------------------------------------------------
## targets


rule all:
    input:
        expand("stages/{sample}.all.done", sample=SAMPLES),


## --------------------------------------------------------------------------------
## rules


rule db_preload:
    output:
        "stages/db.preload.done",
    threads: 64
    benchmark:
        "benchmarks/db.preload.txt"
    shell:
        """          
        krakenuniq --db {DB} --preload --threads {threads} 
        touch {output}
        """


rule classify:
    input:
        fq=get_fq,
        preload="stages/db.preload.done",
    output:
        classify="classify/{sample}." + PREFIX + ".krakenuniq.class.tsv.gz",
        rep="report/{sample}." + PREFIX + ".krakenuniq.report.tsv.gz",
    benchmark:
        "benchmarks/{sample}/classify.txt"
    log:
        "logs/{sample}/classify.log",
    params:
        sample="{sample}",
    threads: 64
    shell:
        """     
        (krakenuniq --db {DB} --threads {threads} --only-classified-output --report-file report/{params.sample}.{PREFIX}.krakenuniq.report.tsv {input.fq} | gzip > {output.classify}) 2> {log}
        gzip report/{params.sample}.{PREFIX}.krakenuniq.report.tsv
        """


rule classify_collect:
    input:
        classify="classify/{sample}." + PREFIX + ".krakenuniq.class.tsv.gz",
        rep="report/{sample}." + PREFIX + ".krakenuniq.report.tsv.gz",
    output:
        stamp="stages/{sample}.classify.done",
    shell:
        """
        touch {output}
        """
        

rule get_genera_sample:
    input:
        reportfile="report/{sample}." + PREFIX + ".krakenuniq.report.tsv.gz",
        genfile=GENFILE,
    output:
        "taxlists/{sample}/genera.txt",
    wildcard_constraints:
        genus="\d+",
    params:
        sample="{sample}",
    shell:
        """
        (gzip -cd {input.reportfile} | mawk 'FNR==NR {{ a[$1]; next }} ($7 in a && $4 >= {KMERS}){{print $7"\\t"$9}}' <(cut -f1 {input.genfile}) -) > {output}
        """


checkpoint get_assemblies_sample:
    input:
        genfile="taxlists/{sample}/genera.txt",
        reportfile="report/{sample}." + PREFIX + ".krakenuniq.report.tsv.gz",
        tax_ids=DB + "/taxLists/genus.taxIds.tsv.gz",
        seq_info=DB + "/library.seqInfo.tsv",
    output:
        directory("taxlists/{sample}/assemblies/"),
        tax_ids="taxlists/{sample}/assemblies/taxIds.txt",
        assembly_ids="taxlists/{sample}/assemblies/assemblyIds.txt",
        assembly_counts="taxlists/{sample}/assemblies/assemblyCounts.txt",
    params:
        sample="{sample}",
    shell:
        """
        gzip -cd {input.tax_ids} | mawk -F'\\t' 'FNR==NR {{ a[$1]; next }} ($2 in a)' <(cut -f1 {input.genfile}) - > {output.tax_ids}
        cat {input.seq_info} | mawk -F'\\t' 'FNR==NR {{ a[$4]; next }} ($2 in a){{print $1"\\t"$3"\\t"$5"\\t"$6}}' {output.tax_ids} - | sort > {output.assembly_ids}
        gzip -cd {input.reportfile} | mawk 'FNR==NR {{ a[$4]; next }} ($7 in a && $8 == "assembly"){{print $9"\\t"$4}}' {output.tax_ids} - | sort > taxlists/{params.sample}/assemblies/assemblyCountsObs.txt
        grep -vFwf <(cut -f1 taxlists/{params.sample}/assemblies/assemblyCountsObs.txt) taxlists/{params.sample}/assemblies/assemblyIds.txt  | mawk '{{print $1"\\t0"}}' > taxlists/{params.sample}/assemblies/assemblyCounts0.txt
        cat taxlists/{params.sample}/assemblies/assemblyCountsObs.txt taxlists/{params.sample}/assemblies/assemblyCounts0.txt | sort > {output.assembly_counts}
        join {output.assembly_ids} {output.assembly_counts} | sort -k3,3 -k2,2 -rnk5,5 | tr ' ' '\\t' | datamash -g2 first 1 first 3 | awk '{{print >"taxlists/{params.sample}/assemblies/"$3"."$2".id"}}'
        cat taxlists/{params.sample}/assemblies/*.id | cut -f3 | sort | uniq | awk '{{print >"taxlists/{params.sample}/assemblies/"$1".genus"}}'
        """


rule get_read_ids_targets:
    input:
        classfile="classify/{sample}." + PREFIX + ".krakenuniq.class.tsv.gz",
        tax_ids="taxlists/{sample}/assemblies/taxIds.txt",
        genfile="taxlists/{sample}/genera.txt",
    output:
        ids=temp(TMP_DIR + "/{sample}/targets.taxIds.txt"),
        classfile=temp(TMP_DIR + "/{sample}/targets.class.tsv.gz"),
        read_ids=temp(TMP_DIR + "/{sample}/targets.readIds.tsv.gz"),
    benchmark:
        "benchmarks/{sample}/get_read_ids_targets.txt"
    threads: 16
    shell:
        """
        cat {input.tax_ids} | awk '{{print $4"\\t"$2}}' > {output.ids}
        cat {input.genfile} | awk '{{print $1"\\t"$1}}' >> {output.ids}
        gzip -cd {input.classfile} | mawk 'FNR==NR {{ a[$1]; next }} ($3 in a){{print $3"\\t"$2}}' <(cut -f1 {output.ids}) - | sort -k1,1 -S 128G --parallel={threads} | gzip > {output.classfile} 
        join <(gzip -cd {output.classfile}) <(sort -k1,1 {output.ids}) | awk '{{print $2"\\t"$3}}' | gzip > {output.read_ids}
        """


rule get_reads_genera_all:
    input:
        fq=get_fq,
        read_ids=TMP_DIR + "/{sample}/targets.readIds.tsv.gz",
    output:
        fq=temp(TMP_DIR + "/{sample}/classified.fq.gz"),
    benchmark:
        "benchmarks/{sample}/get_reads_genera_all.txt"
    priority: 10
    wildcard_constraints:
        genus="\d+",
    shell:
        """
        seqtk subseq <(gzip -cdf {input.fq}) <(gzip -cd {input.read_ids} | cut -f1) | gzip > {output.fq} 
        """


rule get_read_ids_genus:
    input:
        read_ids=TMP_DIR + "/{sample}/targets.readIds.tsv.gz",
    output:
        read_ids=temp(TMP_DIR + "/{sample}/{genus}.readIds.txt.gz"),
    benchmark:
        "benchmarks/{sample}/{genus}.get_read_ids_genus.txt"
    wildcard_constraints:
        genus="\d+",
    shell:
        """
        gzip -cd {input.read_ids} | awk '$2 == {wildcards.genus}{{print $1}}' | gzip > {output.read_ids}
        """


rule get_reads_genus:
    input:
        fq=TMP_DIR + "/{sample}/classified.fq.gz",
        read_ids=TMP_DIR + "/{sample}/{genus}.readIds.txt.gz",
    output:
        fq="fq/{sample}/{genus}.fq.gz",
    benchmark:
        "benchmarks/{sample}/{genus}.get_reads_genus.txt"
    wildcard_constraints:
        genus="\d+",
    shell:
        """
        seqtk subseq {input.fq} {input.read_ids} | seqkit rmdup | gzip > {output.fq} 
        """


rule map_minimap2:
    input:
        fq="fq/{sample}/{genus}.fq.gz",
        mmi=get_mmi,
    output:
        bam=temp(TMP_DIR + "/{sample}/{genus}.{assembly}.init.bam"),
    log:
        "logs/{sample}/{genus}.{assembly}.map_minimap2.log",
    benchmark:
        "benchmarks/{sample}/{genus}.{assembly}.map_minimap2.txt"
    wildcard_constraints:
        genus="\d+",
    threads: 4
    params:
        rg="@RG\\tID:{genus}.{sample}\\tSM:{sample}",
    shell:
        """
        (minimap2 -R '{params.rg}' -t {threads} {MM2_PARAM} -a {input.mmi} {input.fq} | samtools view -F4 -bh > {output.bam}) 2> {log}
        """


rule sort_bam:
    input:
        bam=TMP_DIR + "/{sample}/{genus}.{assembly}.init.bam",
        fa=get_fa,
    output:
        bam=temp(TMP_DIR + "/{sample}/{genus}.{assembly}.srt.bam"),
    wildcard_constraints:
        genus="\d+",
    threads: 2
    shell:
        """
        samtools sort -@{threads} {input.bam} | samtools calmd -b - {input.fa} > {output.bam}
        """


rule mark_duplicates:
    input:
        bam=TMP_DIR + "/{sample}/{genus}.{assembly}.srt.bam",
    output:
        bam="bam/{sample}/{genus}.{assembly}." + PREFIX + ".mrkdup.bam",
        metrics="bam/{sample}/{genus}.{assembly}." + PREFIX + ".mrkdup.metrics.txt",
    wildcard_constraints:
        genus="\d+",
    benchmark:
        "benchmarks/{sample}/{genus}.{assembly}.mark_duplicates.txt"
    log:
        "logs/{sample}/{genus}.{assembly}.mark_duplicates.log",
    threads:
        4
    shell:
        """
        picard -Xmx50g MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} --TMP_DIR {TMP_DIR} 2> {log}
        """


rule index_bam:
    input:
        bam="bam/{sample}/{genus}.{assembly}." + PREFIX + ".mrkdup.bam",
    wildcard_constraints:
        genus="\d+",
    output:
        bai="bam/{sample}/{genus}.{assembly}." + PREFIX + ".mrkdup.bam.bai",
    shell:
        """
        samtools index {input.bam} 
        """


rule filter_bam:
    input:
        bam="bam/{sample}/{genus}.{assembly}." + PREFIX + ".mrkdup.bam",
    output:
        bam="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.bam",
        bai="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.bam.bai",
    wildcard_constraints:
        genus="\d+",
    shell:
        """
        samtools view -q{MQ} -F 0x400 -bh {input.bam} > {output.bam}
        samtools index {output.bam}
        """


rule get_coverage_bam:
    input:
        bam="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.bam",
        bai="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.bam.bai",
    output:
        gc="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.genomecov",
        cov="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.coverage.txt",
    wildcard_constraints:
        genus="\d+",
    shell:
        """
        touch {output.gc}
        (bedtools genomecov -ibam {input.bam} > {output.gc}) || true
        cat {output.gc} | mawk '{{print $1"\\t"$2*$3"\\t"$4}}' | datamash -g1 sum 2 first 3 | mawk '{{print $1,$2/$3}}' > {output.cov}
        """


def aggregate_bams(wildcards):
    checkpoint_output = checkpoints.get_assemblies_sample.get(**wildcards).output[0]
    files = expand(
        "bam/{sample}/{unit}." + PREFIX + ".filter.{ext}",
        sample=wildcards.sample,
        unit=glob_wildcards(os.path.join(checkpoint_output, "{unit}.id")).unit,
        ext=["bam", "bam.bai", "genomecov", "coverage.txt"],
    )
    return files


rule map_collect:
    input:
        aggregate_bams,
    output:
        stamp="stages/{sample}.map.done",
    shell:
        """
        touch {output}
        """


rule sort_bam_metaDMG:
    input:
        bam="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.bam",
        bai="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.bam.bai",
    output:
        bam=TMP_DIR + "/{sample}/{genus}.{assembly}.filter.srt_name.bam",
    wildcard_constraints:
        genus="\d+",
    shell:
        """
        samtools sort -n {input.bam} > {output.bam}
        """

        
rule get_damage_global:
    input:
        bam=TMP_DIR + "/{sample}/{genus}.{assembly}.filter.srt_name.bam",
    output:
        dam=(
            "metadamage/{sample}/{genus}.{assembly}."
            + PREFIX
            + ".filter.damage_global.txt"
        ),
        bdam="metadamage/{sample}/{genus}.{assembly}." + PREFIX + ".filter.global.bdamage.gz",
    benchmark:
        "benchmarks/{sample}/{genus}.{assembly}.get_damage.txt"
    log:
        "logs/{sample}/{genus}.{assembly}.get_damage.log",
    wildcard_constraints:
        genus="\d+",
    params:
        px="metadamage/{sample}/{genus}.{assembly}." + PREFIX + ".filter.global"
    
    shell:
        """
        ({METADMG_CPP} getdamage -p 25 -r 0 {input.bam} -o {params.px}) 2> {log} 1> {output.dam}
        """


rule get_damage_local:
    input:
        bam=TMP_DIR + "/{sample}/{genus}.{assembly}.filter.srt_name.bam",
    output:
        dam=(
            "metadamage/{sample}/{genus}.{assembly}."
            + PREFIX
            + ".filter.damage_local.txt"
        ),
        bdam="metadamage/{sample}/{genus}.{assembly}." + PREFIX + ".filter.local.bdamage.gz",
    benchmark:
        "benchmarks/{sample}/{genus}.{assembly}.get_damage.txt"
    log:
        "logs/{sample}/{genus}.{assembly}.get_damage.log",
    wildcard_constraints:
        genus="\d+",
    params:
        px="metadamage/{sample}/{genus}.{assembly}." + PREFIX + ".filter.local"
    
    shell:
        """
        ({METADMG_CPP} getdamage -p 25 -r 1 {input.bam} -o {params.px}) 2> {log} 1> {output.dam}
        """


def aggregate_damage(wildcards):
    checkpoint_output = checkpoints.get_assemblies_sample.get(**wildcards).output[0]
    files = expand(
        "metadamage/{sample}/{unit}." + PREFIX + ".filter.{ext}",
        sample=wildcards.sample,
        unit=glob_wildcards(os.path.join(checkpoint_output, "{unit}.id")).unit,
        ext=["global.bdamage.gz", "damage_global.txt", "local.bdamage.gz", "damage_local.txt"],
    )
    return files


rule damage_collect:
    input:
        aggregate_damage,
    output:
        stamp="stages/{sample}.damage.done",
    shell:
        """
        touch {output}
        """


rule summary_sample:
    input:
        "stages/{sample}.map.done",
        "stages/{sample}.damage.done",
    output:
        tsv="tables/{sample}/" + PREFIX + ".summary.tsv.gz",
        stamp="stages/{sample}.summary.done",
    params:
        sample="{sample}",
    threads: 48
    shell:
        """
        Rscript src/summary.R {PREFIX} {DB} {params.sample} {threads} {METADMG_CPP}
        touch {output.stamp}
        """


def aggregate_bams_genus(wildcards):
    checkpoint_output = checkpoints.get_assemblies_sample.get(**wildcards).output[0]
    units = glob_wildcards(os.path.join(checkpoint_output, "{unit}.id")).unit
    g = [u.split(".")[0] for u in units]
    a = [u.split(".", 1)[1] for u in units]
    f = [a[i] for i in range(len(a)) if g[i] == wildcards.genus]
    files = expand(
        "bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.{ext}",
        sample=wildcards.sample,
        genus=wildcards.genus,
        assembly=f,
        ext=["bam"],
    )
    return files


rule plot_edit_dist:
    input:
        aggregate_bams_genus,
    output:
        pdf="plots/{sample}/{genus}." + PREFIX + ".editDist.pdf",
        tsv="tables/{sample}/{genus}." + PREFIX + ".editDist.tsv",
    wildcard_constraints:
        genus="\d+",
    shell:
        """
        Rscript src/plotEditDist.R {PREFIX} {DB} {output.pdf} {output.tsv} {input}
        """


def aggregate_damage_genus(wildcards):
    checkpoint_output = checkpoints.get_assemblies_sample.get(**wildcards).output[0]
    units = glob_wildcards(os.path.join(checkpoint_output, "{unit}.id")).unit
    g = [u.split(".")[0] for u in units]
    a = [u.split(".", 1)[1] for u in units]
    f = [a[i] for i in range(len(a)) if g[i] == wildcards.genus]
    files = expand(
        "metadamage/{sample}/{genus}.{assembly}." + PREFIX + ".filter.{ext}",
        sample=wildcards.sample,
        genus=wildcards.genus,
        assembly=f,
        ext=["global.bdamage.gz"],
    )
    return files


rule plot_damage:
    input:
        aggregate_damage_genus,
    output:
        pdf="plots/{sample}/{genus}." + PREFIX + ".damage.pdf",
    wildcard_constraints:
        genus="\d+",
    shell:
        """
        Rscript src/plotDamage.R {PREFIX} {DB} {output.pdf} {METADMG_CPP} {input}
        """


def aggregate_plots(wildcards):
    checkpoint_output = checkpoints.get_assemblies_sample.get(**wildcards).output[0]
    files = expand(
        "plots/{{sample}}/{genus}." + PREFIX + ".{ext}",
        genus=glob_wildcards(os.path.join(checkpoint_output, "{genus}.genus")).genus,
        ext=["damage.pdf", "editDist.pdf"],
    )
    return files


rule plots_collect:
    input:
        aggregate_plots,
    output:
        stamp="stages/{sample}.plots.done",
    shell:
        """
        touch {output}
        """


rule all_collect:
    input:
        "stages/{sample}.classify.done",
        "stages/{sample}.map.done",
        "stages/{sample}.plots.done",
        "stages/{sample}.summary.done",
    output:
        stamp="stages/{sample}.all.done",
    shell:
        """
        touch {output}
        """
