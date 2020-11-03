
import pandas as pd
import glob


## --------------------------------------------------------------------------------
## global parameters from config file

PREFIX=config["prefix"] ## output prefix
DB=config["db"] ## database root dir
GENFILE=config["genera"] ## target genera file
KMERS=config["kmers"] ## minimum kmers for target assemblies
MQ=config["MQ"] ## MQ cutoff for mapped BAMs


## --------------------------------------------------------------------------------
## helpers

unit_df = pd.read_table(config["units"], comment="#").set_index(["sampleId"], drop=False)
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
        expand("stages/{sample}.all.done", sample = SAMPLES)

          
## --------------------------------------------------------------------------------
## rules

rule db_preload:
    output:
        "stages/db.preload.done"
    threads:
        48
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
        preload="stages/db.preload.done"
    output:
        classify="classify/{sample}." + PREFIX + ".krakenuniq.class.tsv.gz",
        rep="report/{sample}." + PREFIX + ".krakenuniq.report.tsv.gz",
    benchmark:
        "benchmarks/{sample}/classify.txt"
    log:
        "logs/{sample}/classify.log"
    params:
        sample="{sample}"
    threads:
        48
    shell:
        """     
        (krakenuniq --db {DB} --threads {threads} --only-classified-output --report-file report/{params.sample}.{PREFIX}.krakenuniq.report.tsv {input.fq} | gzip > {output.classify}) 2> {log}
        gzip report/{params.sample}.{PREFIX}.krakenuniq.report.tsv
        """

rule get_read_ids_classified_nonhuman:
    input:
        classfile="classify/{sample}." + PREFIX + ".krakenuniq.class.tsv.gz",
        tax_ids_human=DB + "/taxLists/33208.kingdom.taxIds.txt"
    output:
        read_ids="tmp/{sample}/classified.readIds.txt.gz"
    benchmark:
        "benchmarks/{sample}/get_read_ids_classified.txt"
    shell:
        """
        gzip -cd {input.classfile} | awk 'FNR==NR {{ a[$1]; next }} (!($3 in a) && $1 == "C"){{print $2}}' {input.tax_ids_human} - | gzip > {output.read_ids} 
        """

rule get_reads_classified_nonhuman:
    input:
        fq=get_fq,
        read_ids="tmp/{sample}/classified.readIds.txt.gz"
    output:
        fq="fq/{sample}/classified.fq.gz"
    benchmark:
        "benchmarks/{sample}/get_reads_classified_nonhuman.txt"
    threads:
        4
    shell:
        """
        seqkit grep -f {input.read_ids} {input.fq} | pigz -p {threads} > {output.fq} 
        """

rule get_genera_sample:
    input:
        reportfile="report/{sample}." + PREFIX + ".krakenuniq.report.tsv.gz",
        genfile=GENFILE
    output:
        "tmp/{sample}/genera.txt"
    wildcard_constraints:
        genus="\d+"
    params:
        sample="{sample}",
    shell:
        """
        (gzip -cd {input.reportfile} | awk 'FNR==NR {{ a[$1]; next }} ($7 in a && $4 >= {KMERS}){{print $7"\t"$9}}' <(cut -f1 {input.genfile}) -) > {output}
        """
        
checkpoint get_assemblies_sample:
    input:
        genfile="tmp/{sample}/genera.txt",
        reportfile="report/{sample}." + PREFIX + ".krakenuniq.report.tsv.gz",
        tax_ids=DB + "/taxLists/genus.taxIds.tsv.gz",
        seq_info=DB + "/library.seqInfo.tsv",
    output:
        directory("tmp/{sample}/assemblies/"),
        tax_ids="tmp/{sample}/assemblies/taxIds.txt",
        assembly_ids="tmp/{sample}/assemblies/assemblyIds.txt",
        assembly_counts="tmp/{sample}/assemblies/assemblyCounts.txt"
    params:
        sample="{sample}",
    shell:
        """
        zcat {input.tax_ids} | awk -F'\\t' 'FNR==NR {{ a[$1]; next }} ($2 in a){{print $4}}' <(cut -f1 {input.genfile}) - > {output.tax_ids}
        cat {input.seq_info} | awk -F'\\t' 'FNR==NR {{ a[$1]; next }} ($2 in a){{print $1"\\t"$3"\\t"$5}}' {output.tax_ids} - | sort > tmp/{params.sample}/assemblies/assemblyIds.txt
        gzip -cd {input.reportfile} | awk 'FNR==NR {{ a[$1]; next }} ($7 in a && $8 == "assembly"){{print $9"\\t"$4}}' {output.tax_ids} - | sort > tmp/{params.sample}/assemblies/assemblyCounts.txt
        if [ -s tmp/{params.sample}/assemblies/assemblyCounts.txt ]
        then
           join tmp/{params.sample}/assemblies/assemblyIds.txt tmp/{params.sample}/assemblies/assemblyCounts.txt | sort -k2,2 -rnk3,3 | tr ' ' '\\t' | datamash -g2 first 1 first 3 | awk '{{print >"tmp/{params.sample}/assemblies/"$3"."$2".id"}}'
        else
           cat tmp/{params.sample}/assemblies/assemblyIds.txt | awk '{{print $1":"$3"\\t"$2}}' | datamash -g2 rand 1 | tr ':' '\\t' | awk '{{print >"tmp/{params.sample}/assemblies/"$3"."$2".id"}}'
        fi
        cat tmp/{params.sample}/assemblies/*.id | cut -f2 | sort | uniq | awk '{{print >"tmp/{params.sample}/assemblies/"$1".assembly"}}'
        cat tmp/{params.sample}/assemblies/*.id | cut -f3 | sort | uniq | awk '{{print >"tmp/{params.sample}/assemblies/"$1".genus"}}'
        """
        
rule get_read_ids_genus:
    input:
        classfile="classify/{sample}." + PREFIX + ".krakenuniq.class.tsv.gz",
        genera="tmp/{sample}/assemblies/{genus}.genus",
        tax_ids=DB + "/taxLists/{genus}.genus.taxIds.txt"
    output:
        read_ids=temp("tmp/{sample}/{genus}.readIds.txt.gz")
    benchmark:
        "benchmarks/{sample}/{genus}.get_read_ids_genus.txt"
    shell:
        """
        gzip -cd {input.classfile} | awk 'FNR==NR {{ a[$1]; next }} ($3 in a){{print $2}}' {input.tax_ids} - | gzip > {output.read_ids}
        """

rule get_reads_genus:
    input:
        fq="fq/{sample}/classified.fq.gz",
        read_ids="tmp/{sample}/{genus}.readIds.txt.gz"
    output:
        fq="fq/{sample}/{genus}.fq.gz"
    benchmark:
        "benchmarks/{sample}/{genus}.get_reads_genus.txt"
    shell:
        """
        seqkit grep -f {input.read_ids} {input.fq} | gzip > {output.fq} 
        """

# def aggregate_fqs(wildcards):
#     checkpoint_output = checkpoints.get_assemblies_sample.get(**wildcards).output[0]    
#     files = expand("fq/{{sample}}/{genus}.fq.gz", genus = glob_wildcards(os.path.join(checkpoint_output, "{genus}.genus")).genus)
#     return files 

# rule fqs_collect:
#     input:
#         aggregate_fqs
#     output:
#         stamp = "stages/{sample}.fqs.done"
#     shell:
#         """
#         touch {output}
#         """

rule map_minimap2:
    input:
        fq="fq/{sample}/{genus}.fq.gz",
        mmi=get_mmi,
    output:
        bam=temp("tmp/{sample}/{genus}.{assembly}.init.bam"),
    log:
        "logs/{sample}/{genus}.{assembly}.map_minimap2.log"
    benchmark:
        "benchmarks/{sample}/{genus}.{assembly}.map_minimap2.txt"
    wildcard_constraints:
        genus="\d+"
    threads:
        4
    params:
        rg="@RG\\tID:{genus}.{sample}\\tSM:{sample}"
    shell:
        """
        (minimap2 -R '{params.rg}' -t {threads} -x sr -a {input.mmi} {input.fq} | samtools view -F4 -bh > {output.bam}) 2> {log}
        """

rule sort_bam:
    input:
        bam="tmp/{sample}/{genus}.{assembly}.init.bam"
    output:
        bam=temp("tmp/{sample}/{genus}.{assembly}.srt.bam"),
    wildcard_constraints:
        genus="\d+"
    threads:
        2
    shell:
        """
        samtools sort -@{threads} {input.bam} > {output.bam}
        """
        
rule mark_duplicates:
    input:
        bam="tmp/{sample}/{genus}.{assembly}.srt.bam",
    output:
        bam="bam/{sample}/{genus}.{assembly}." + PREFIX + ".mrkdup.bam",
        metrics="bam/{sample}/{genus}.{assembly}." + PREFIX + ".mrkdup.metrics.txt"
    wildcard_constraints:
        genus="\d+"
    benchmark:
        "benchmarks/{sample}/{genus}.{assembly}.mark_duplicates.txt"
    log:
        "logs/{sample}/{genus}.{assembly}.mark_duplicates.log"
    shell:
        """
        java -Xmx16g -jar /willerslev/software/picard/picard.jar MarkDuplicates I={input} O={output.bam} M={output.metrics} 2> {log}
        """

rule index_bam:
    input:
        bam="bam/{sample}/{genus}.{assembly}." + PREFIX + ".mrkdup.bam",
    wildcard_constraints:
        genus="\d+"
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
        bai="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.bam.bai"
    wildcard_constraints:
        genus="\d+"
    shell:
        """
        samtools view -q{MQ} -F 0x400 -bh {input.bam} > {output.bam}
        samtools index {output.bam}
        """

rule get_coverage_bam:
    input:
        bam="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.bam",
        bai="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.bam.bai"
    output:
        gc="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.genomecov",
        cov="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.coverage.txt"
    wildcard_constraints:
        genus="\d+"
    shell:
        """
        touch {output.gc}
        (bedtools genomecov -ibam {input.bam} > {output.gc}) || true
        cat {output.gc} | awk '{{print $1"\\t"$2*$3"\\t"$4}}' | datamash -g1 sum 2 first 3 | awk '{{print $1,$2/$3}}' > {output.cov}
        """

def aggregate_bams(wildcards):
    checkpoint_output = checkpoints.get_assemblies_sample.get(**wildcards).output[0]
    files = expand("bam/{{sample}}/{unit}." + PREFIX + ".filter.{ext}", unit = glob_wildcards(os.path.join(checkpoint_output, "{unit}.id")).unit, ext = ['bam', 'bam.bai', 'genomecov', 'coverage.txt'])
    return files 


rule map_collect:
    input:
        aggregate_bams
    output:
        stamp="stages/{sample}.map.done"
    shell:
        """
        touch {output}
        """
          
rule get_damage:
    input:
        bam="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.bam",
        bai="bam/{sample}/{genus}.{assembly}." + PREFIX + ".filter.bam.bai",
        fa=get_fa
    output:
        mis="mapdamage/{sample}/{genus}.{assembly}." + PREFIX + ".filter/misincorporation.txt.gz",
    benchmark:
        "benchmarks/{sample}/{genus}.{assembly}.get_damage.txt"
    params:
        sample="{sample}",
        genus="{genus}",
        assembly="{assembly}",
    wildcard_constraints:
        genus="\d+"
    shell:
        """
        mapDamage -i {input.bam} -r {input.fa} -d mapdamage/{params.sample}/{params.genus}.{params.assembly}.{PREFIX}.filter --no-stats
        gzip mapdamage/{params.sample}/{params.genus}.{params.assembly}.{PREFIX}.filter/misincorporation.txt
        """

def aggregate_damage(wildcards):
    checkpoint_output = checkpoints.get_assemblies_sample.get(**wildcards).output[0]
    files = expand("mapdamage/{{sample}}/{unit}." + PREFIX + ".filter/misincorporation.txt.gz", unit = glob_wildcards(os.path.join(checkpoint_output, "{unit}.id")).unit)
    return files 

rule damage_collect:
    input:
        aggregate_damage
    output:
        stamp="stages/{sample}.damage.done"
    shell:
        """
        touch {output}
        """

rule summary_sample:
    input:
        "stages/{sample}.map.done",
        "stages/{sample}.damage.done",
    output:
        tsv="tables/{sample}/" + PREFIX + ".summary.tsv",
        stamp = "stages/{sample}.summary.done"
    params:
        sample="{sample}"
    threads:
        12
    shell:
        """
        Rscript src/summary.R {PREFIX} {DB} {params.sample} {threads}
        touch {output.stamp}
        """
        
rule plot_edit_dist:
    input:
        stamp="stages/{sample}.map.done",
        files=lambda wildcards: glob.glob("bam/{sample}/{genus}.*filter.bam".format(sample=wildcards.sample, genus=wildcards.genus))
    output:
        pdf="plots/{sample}/{genus}." + PREFIX + ".editDist.pdf",
        tsv="tables/{sample}/{genus}." + PREFIX + ".editDist.tsv",
    wildcard_constraints:
        genus="\d+"
    shell:
        """
        Rscript src/plotEditDist.R {DB} {output.pdf} {output.tsv} {input.files}
        """

rule plot_damage:
    input:
        stamp="stages/{sample}.damage.done",
        files=lambda wildcards: glob.glob("mapdamage/{sample}/{genus}.*filter/misincorporation.txt.gz".format(sample=wildcards.sample, genus=wildcards.genus))
    output:
        pdf="plots/{sample}/{genus}." + PREFIX + ".damage.pdf",
        tsv="tables/{sample}/{genus}." + PREFIX + ".damage.tsv.gz",
    wildcard_constraints:
        genus="\d+"
    shell:
        """
        Rscript src/plotDamage.R {DB} {output.pdf} {output.tsv} {input.files}
        """

def aggregate_plots(wildcards):
    checkpoint_output = checkpoints.get_assemblies_sample.get(**wildcards).output[0]    
    files = expand("plots/{{sample}}/{genus}." + PREFIX + ".{ext}", genus = glob_wildcards(os.path.join(checkpoint_output, "{genus}.genus")).genus, ext = ['damage.pdf', 'editDist.pdf'])
    return files 

rule plots_collect:
    input:
        aggregate_plots
    output:
        stamp="stages/{sample}.plots.done"
    shell:
        """
        touch {output}
        """  

rule all_collect:
    input:
        "stages/{sample}.plots.done",
        "stages/{sample}.summary.done"
    output:
        stamp="stages/{sample}.all.done"
    shell:
        """
        touch {output}
        """  
