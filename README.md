# Prok-SNPTree

[![DOI](https://zenodo.org/badge/423667401.svg)](https://zenodo.org/badge/latestdoi/423667401)

Prok-SNPTree is a pipeline to generate phylogenetic tree using core substitutions between reference genome and whole genomes Illumina sequencing.  
Prok-SNPTree was originally designed and optimized for small prokaryotic genomes, but can perform reasonably well with larger genomes up to 100 Mb.

Prok-SNPtree comes with a helper file for providing required input info. The sbatch has preset information for running through SLURM scheduler, but should perform without SLURM as long as the following executables are available in the ENV.

| Function | Tools/scripts |
| --- | --- |
| Parallelized sample processing | **GNU Parallel** ⇨ [Source](https://gnu.askapache.com/parallel/) · [Website](https://www.gnu.org/software/parallel/) |
| Quality check for raw reads | **FᴀsᴛQC** ⇨ [Source](https://github.com/s-andrews/FastQC) · [Website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) <br /> **MᴜʟᴛɪQC** ⇨ [Source](https://github.com/ewels/MultiQC) · [Reference](http://dx.doi.org/10.1093/bioinformatics/btw354) · [Website](https://multiqc.info/) |
| Adapter identification and trimming | **ᴄᴜᴛᴀᴅᴀᴘᴛ** ⇨ [Source](https://github.com/marcelm/cutadapt/) · [Reference](http://dx.doi.org/10.14806/ej.17.1.200) <br /> **ᴛʀɪᴍ_ɢᴀʟᴏʀᴇ** ⇨ [Source](https://github.com/FelixKrueger/TrimGalore) · [Reference](https://doi.org/10.5281/zenodo.5127899) · [Website](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) |
| Genome indexing and read alignment | **ʙᴡᴀ** ⇨ [Source](https://github.com/lh3/bwa) · [Reference](https://doi.org/10.1093/bioinformatics/btp324) · [Reference](https://doi.org/10.48550/arXiv.1303.3997) · [Website](https://bio-bwa.Sourceforge.net/) |
| Binary conversion and sorting | **Sᴀᴍᴛᴏᴏʟs** ⇨ [Source](https://github.com/samtools/samtools) · [Reference](https://doi.org/10.1093/bioinformatics/btp352) · [Reference](https://doi.org/10.1093/gigascience/giab008) · [Website](http://www.htslib.org/) |
| Variant calling and selection | **GATK** ⇨ [Source](https://github.com/broadinstitute/gatk) · [Reference](http://dx.doi.org/10.1038/ng.806) · [Website](https://gatk.broadinstitute.org/) |
| Phylogenetic tree | **RAxML** ⇨ [Source](https://github.com/stamatak/standard-RAxML) · [Reference](https://doi.org/10.1093/bioinformatics/btu033) · [Website](https://cme.h-its.org/exelixis/web/software/raxml/) |
| Pairwise SNP count | **FᴀsᴛᴀTᴏSNPCᴏᴜɴᴛ.sʜ** ⇨ [Source](https://gist.github.com/rknx/3d3ad3b93ad963be84d7f2840486e07f) |

### Input files
- All paired gzipped fastq can be placed in working directory or in a subdirectory named `fastq`. The pipeline does not support singletons currently.
- Reference genome should be in `refs` subdirectory, and named `genome.fna`. Symbolic links are accepted. Alternatively, the script can download it automatically (with wget) if a direct link is provided (see arguments below).
- Reference annotation is not currently. For futureproofing, it may be provided inside `refs` subdirectory as `genes.gtf`. See alternative methods in arguments.

### Arguments

The script accepts the following arguments, which are supplied from the helper sbatch file.

1. **refgenome** (Reference genome)  
Direct link (url) to reference genome (gzipped). This for convenience, and the intended goal is to be able to download reference genome from NCBI etc. Keep its value empty to `.` if it will be supplied maunally (see input files above).

2. **refannotation** (Annotation file for reference genome)  
Direct link (url) to reference genome (gzipped) as `refgenome`. This option is here for future function, and may be set to empty now.

3. **minDP**, **minQD**, **maxRDP**, and **minADP** (VCF filtration parameters)

    - minDP: Minimum sequencing depth (Positions that fail are considered absent)
    - minQD: Minimum depth-normalized quality (SNPs that fail are ignored)
    - maxRDP: Maximum reference allele depth for SNPs to be accepted as real.
    - minADP: Minimum alternate allele depth for SNPs to be accepted as real.

4. **ncpu** and **nmem** (Parallelization parameters)  
Number of CPUs and memory (overall) to use. If number of CPU is less than 8, only one sample is processed at a time with the number of CPU available. Otherwise, 8 CPUs are used per sample, and the samples are parallelized based on number of CPUs.

5. **nboot** (Bootstrap)  
Number of bootstraps to be used while preparing the phylogenetic tree with RAxML.

### Optimization
The pipeline is written so as to enable resuming or rerunning. The main output files are name systematically and are used as checkpoints.  
Some examples:
- Trimming: non-empty `\<sample>.og` files in **fastq** subdirectory.
- Alignment: non-empty `.bam` file in **align/\<sample>** subdirectory.
- Variant calling: non-empty `.vcf` file in **variants/\<sample>** subdirectory.

> If some samples fail, just rerun the pipeline after amking changes with the inputs. Completed samples with valid outputs are not processed again.  

> If the pipeline exits in middle of operation, just rerun, and the pipeline will pick up from last complete operation.

> If you add new samples, just rerun the pipeline to process it.

### Citation
A peer-reviewed paper is pending publication. Please cite the zenodo record at the moment as follows:  
> Sharma A. 2022. rknx/Prok-SNPTree: Phylogenetic Tree from Next-gen Sequencing of Prokaryotes (v0.1b). Zenodo. https://doi.org/10.5281/zenodo.7444860

