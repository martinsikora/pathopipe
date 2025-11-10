# pathopipe
**Ancient Pathogen Screening Workflow**

The pipeline is described in detail in Sikora et al., 2025, **"The spatiotemporal distribution of human pathogens in ancient Eurasia"** [https://www.nature.com/articles/s41586-025-09192-8](https://www.nature.com/articles/s41586-025-09192-8)

## Table of Contents
- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
- [Output Files](#output-files)
- [License](#license)

## About
`pathopipe` is a Snakemake workflow designed to identify and classify microbial DNA within ancient shotgun sequencing data. The workflow was developed to detect pathogens in ancient human data, but it can be applied across a wide range of microbial and eukaryotic targets, using sequencing data from various sources such as animal remains and ancient environmental samples. 

The pipeline orchestrates multiple steps—quality control, mapping, damage profiling, taxonomic classification, summary reporting—into a unified, reproducible framework.

## Installation
1. Clone the repository:

```
git clone https://github.com/martinsikora/pathopipe.git 
cd pathopipe
```

2. Install required dependencies:
    - krakenuniq (1.0.4)
    - mawk (1.3.4)
    - seqtk (1.3-r106)
    - seqkit (2.3.0)
    - bowtie2 (2.5.2)
    - samtools (1.17)
    - picard (2.27.5)
    - bedtools (2.30.0)
    - datamash (1.5)
    - metaDMG (0.2-41-gc867207)
    - snakemake (7.20.0)
    - gargammel (1.1.4)
    - R (4.2.2)
    - R package fastTopics (0.6-142)
    - R package Rbeast (0.9.7)
    - R package tidyverse (1.3.2)
    - R package inlabru (2.8.0)
    - R package Rsamtools (2.14.0)

If conda is pre-installed on your system, this can be conveniently installed in a conda (or mamba or similar) environment using the following:

   ```
    conda create -n pathopipe -c conda-forge -c bioconda -c defaults \
      snakemake \
      krakenuniq=1.0.4 \
      mawk=1.3.4 \
      seqtk=1.3 \
      seqkit=2.3.0 \
      bowtie2=2.5.2 \
      samtools=1.17 \
      picard=2.27.5 \
      bedtools=2.30.0 \
      datamash=1.5 \
      snakemake=7.20.0 \
      gargammel=1.1.4 \
      metaDMG \
      r-base=4.2.2 \
      r-fasttopics=0.6_142 \
      r-rbeast=0.9.7 \
      r-tidyverse=1.3.2 \
      r-inlabru=2.8.0 \
      bioconductor-rsamtools=2.14.0
   ```
To load the conda environment before running the pipeline, use the following:
   ```    
    conda activate pathopipe
   ```

4. Download and unpack reference database:
   ```
   wget https://erda.ku.dk/archives/1d29e091b69cabe43093440eeb396212/diseases/public_supplementary_data_repo/hum_microbe_release_20250428.tar.gz
   tar -xvzf hum_microbe_release_20250428.tar.gz
   ```

## Usage
Create a tab-separated list of sample-IDs and corresponding fastq files with column names `sampleId` and `fq` (see example file `units.tsv`). Edit your config.yml to point to your units file, reference databases (e.g. `./hum_microbe_release_20250428`), and modify any other parameters if needed.

To run the workflow:

```
snakemake --configfile config.yml --cores <N>
```
Replace `<N>` with the number of CPU cores you wish to allocate.

(Optional) To summarise results across all samples analysed after completion of the pathopipe pipeline, run:
```
snakemake -s summarize.Snakefile --configfile config.yml --cores <N>
```

## Output files
For each sample listed in your units file, a summary table will be generated in `tables/<SAMPLE>/<PREFIX>.summary.tsv.gz`. This output table contains summary statistics for all species within the genera detected. Furthermore, for each sample/genus combination edit distance and damage plots will be generated: `plots/<SAMPLE>/<GENUSTAXID>.<PREFIX>.editDist.pdf` and `plots/<SAMPLE>/<GENUSTAXID>.<PREFIX>.damage.pdf`. 

## License

`pathopipe` is released under the MIT License.

Copyright (c) 2025 Martin Sikora

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

