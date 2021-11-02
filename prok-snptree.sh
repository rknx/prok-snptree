#!/bin/sh

## Prepare environment ---------------------------------------------------------

# Timestamp function
__ts() { echo -e [`date "+%m/%d %H:%M:%S"`]: $1; }
export -f __ts

# Get arguments
export refgenome=$1
export refannotation=$2
export ncpu=$3
export nmem=$4
export qd=$5
export dp=$6
export rdp=$7
export adp=$8
export boot=$9

# Workspace
__ts "Creating workspace sub-folders."
mkdir fastq tmp align variants output

# Set tmp directory
export TMPDIR=$(pwd)/tmp/

# Setting memory limits for java
export _JAVA_OPTIONS=-Xmx"$(( $nmem / $(($ncpu/8)) ))"g

# Move fastq files to right place
mv *.fastq.gz fastq/

## Download and index reference genome -----------------------------------------

# Load required modules
module -q reset; module load bwa gatk samtools

__ts "Indexing the reference genome."

# Function to index prokaryotic genome
__index_genome() {

    # Reference genome
    [[ ! -s refs/genome.fna ]] && \
        __ts "Downloading and extracting reference genome." && ( \
            [[ `wget --spider $refgenome 2> /dev/null` -eq 0 ]] && \
                wget -q -O - $refgenome | gzip -d > refs/genome.fna || \
                __ts "Reference genome link is not valid." \
        ) || \
        __ts "Reference genome is already present.\n\t\tRun 'rm -r refs/*' \
            before this script to rebuild reference files."
    
    # Annotation file
    [[ ! -s refs/genes.gtf ]] && \
        __ts "Downloading reference annotation/" && ( \
            [[ `wget --spider $refannot 2> /dev/null` -eq 0 ]] && \
                wget -q -O - $refannot | gzip -d | grep -v "unknown_transcript" |\
                    awk 'FS=OFS="\t" {if ($4>$5) {s=$4; $4=$5; $5=s}; print }' > \
                        refs/genes.gtf || \
                __ts "Reference annotation link is not valid." \
        ) || __ts "Reference annotation is already present. Skipping."

    # Indexing the genome
    [[ ! -s refs/genome.fna.bwt ]] && \
        __ts "Starting bwa indexing." && \
        bwa index refs/genome.fna || \
        __ts "Reference genome is already indexed. Skipping indexing."

    # Create reference image (this needs .fa or .fasta extension)
    [[ ! -s refs/genome.img ]] && \
        __ts "Creating reference genome image." && \
        ln -s $(pwd)/refs/genome.fna refs/genome.fasta && \
        gatk BwaMemIndexImageCreator -verbosity WARNING \
            -I refs/genome.fasta \
            -O refs/genome.img && \
        unlink refs/genome.fasta || \
        __ts "Reference genome image is already present. Skipping."

    # Generate 2bit file
    [[ ! -s refs/fna22bit || ! -s refs/genome.2bit ]] && \
        __ts "Generating genome 2bit." && \
        wget -q -O - "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit" \
            > refs/fna22bit && \
        chmod +x refs/fna22bit && \
        refs/fna22bit refs/genome.fna refs/genome.2bit || \
        __ts "Genome 2bit is already present. Skipping."

    # Generate bad kmer list
    [[ ! -s refs/bad_kmers.txt ]] && \
        __ts "Generating list of bad kmers." && \
        samtools faidx refs/genome.fna && \
        gatk CreateSequenceDictionary --VERBOSITY WARNING \
            -R refs/genome.fna \
            -O refs/genome.dict && \
        gatk FindBadGenomicKmersSpark -verbosity WARNING \
            -R refs/genome.fna \
            -O refs/bad_kmers.txt || \
        __ts "Bad kmer list is already present. Skipping."

}

# Run the indexing function
__index_genome

## Create custom snpEff annotation database ------------------------------------

# Load required modules
module -q reset; module load snpeff

__ts "Creating snpEff database."

# Function to create snpEff database
__snpeff_dbmake() {
    
    # Make custom snpEff database
    [[ ! -s annot/snpeff/genes.gtf ]] && \
        __ts "Getting reference annotation." && \
        mkdir -p annot/snpeff && \
        ln -s $(pwd)/refs/genes.gtf annot/snpeff/genes.gtf || \
        __ts "snpEff annotation file is already present. Skipping."
    
    # Make custom snpEff database
    [[ ! -s annot/snpeff/sequences.fa ]] && \
        __ts "Getting reference genome." && \
        ln -s $(pwd)/refs/genome.fna annot/snpeff/sequences.fa || \
        __ts "snpEff genome is already present. Skipping."
    
    # Make custom snpEff database
    # Remove codontable for eukaryotes
    [[ ! -s annot/snpeff/snpEffectPredictor.bin ]] && \
        __ts "Preparing database." && \
        snpEff build snpeff \
            -configOption snpeff.genome=org \
            -configOption snpeff.codonTable=Bacterial_and_Plant_Plastid \
            -dataDir $(pwd)/annot || \
        __ts "snpEff database is already present. Skipping."

}

# Run the databse creation function
__snpeff_dbmake

## Merge filepart and lanes ----------------------------------------------------

# Load required modules
module -q reset; module load parallel

__ts "Merging multiple lanes if necessary."

# Standardize the name if necessary, e.g., to add lane infomation:
# find ./ -name "*_R?_00?.fastq.gz" -and -not -name "*_L00?*" -exec \
#     sh -c 'x="{}"; mv "$x" "${x//_R?/_L001_R?}' \;

# Get fastq samples
smlist=( $(
	find fastq/ -name "*R?_00?.fastq.gz" | \
        sed 's/^fastq\///; s/\(_L00.\)\{0,1\}_R._00.\.fastq\.gz$//' | \
        sort -u
) )
__ts "${#smlist[@]} samples found."

__merge_lanes(){

	# Merging lanes if multiple lanes present.
	[[ ! -s fastq/"$1"_R1.fastq.gz && ! -s fastq/"$1"_R2.fastq.gz ]] && \
		__ts "Merging lanes and fileparts for $1." && \
		cat fastq/"$1"_*R1_00?.fastq.gz > fastq/"$1"_R1.fastq.gz && \
		cat fastq/"$1"_*R2_00?.fastq.gz > fastq/"$1"_R2.fastq.gz && \
		find fastq/ -name "$1"_L00?_R?_001.fastq.gz -exec \
			sh -c 'x="{}"; mv "$x" "${x}.bak"' \; || \
        __ts "Lanes and fileparts are already merged for $1. Skipping."

}

# Make function global
export -f __merge_lanes

# Run the merge function
parallel -j $(( $ncpu > 8 ? 8 : $ncpu )) __merge_lanes ::: ${smlist[@]}

__ts "Lane merging completed."

## Check quality ---------------------------------------------------------------

__ts "Checking quality of fastq files."

# Load modules
module -q reset; module load parallel fastqc multiqc

# Get fastq samples
fqlist=( $(
	find fastq/ -name "*_R?.fastq.gz" | \
		sed 's/fastq\///; s/\.fastq\.gz$//'
) )
__ts "${#fqlist[@]} fastq files found."

__check_qual(){

	# Check quality of individual files
	[[ ! -s fastq/"$1"_fastqc.html ]] && \
		__ts "Checking $1." && \
		fastqc -q -t $(( $ncpu > 8 ? $(($ncpu/8)) : 1 )) \
            -o fastq/$(dirname $1) \
            fastq/$1.fastq.gz || \
		__ts "QC result already present. Skipping."

}

# Make function global
export -f __check_qual

# Run the qc function
parallel -j $(( $ncpu > 8 ? 8 : $ncpu )) __check_qual ::: ${fqlist[@]}

# Combine all QC into one file
[[ ! -s fastq/multiqc_report.html ]] && \
    __ts "Generating combines QC report." && \
	multiqc -q -o fastq/ fastq/ || \
	__ts "MultiQC report already exists. Skipping."

# View multiQC output manually or use lynx browser
# The trimming step below should be done manually if necessary.

__ts "QC check completed."

## Trim adapters ---------------------------------------------------------------

__ts "Trimming the fastq files."

# Load modules
module -q reset; module load parallel trim_galore

# Get fastq samples
fqlist=( $(
	find fastq/ -name "*R?.fastq.gz" | \
		sed 's/^fastq\///; s/_R.\.fastq\.gz$//' | \
		sort -u
) )
__ts "${#fqlist[@]} samples found."

# Trimming function wrapper
__trim_adapt(){

	# Check quality of individual files
	[[ ! -s fastq/"$1"_R1.fastq.gz.og && ! -s fastq/"$1"_R2.fastq.gz.og ]] && \
		__ts "Trimming $1." && \
		trim_galore -j $(( $ncpu > 8 ? 8 : $ncpu )) --paired --suppress_warn \
			-o fastq/$(dirname $1) \
			fastq/"$1"_R1.fastq.gz \
			fastq/"$1"_R2.fastq.gz || \
		__ts "$1 seems to be trimmed already. Skipping."

	# Swap trimmed and untrimmed files
	[[ -s fastq/"$1"_R1_val_1.fq.gz && -s fastq/"$1"_R2_val_2.fq.gz ]] && \
		__ts "Backing up original FASTQ files for $1." && \
		mv fastq/"$1"_R1.fastq.gz fastq/"$1"_R1.fastq.gz.og && \
		mv fastq/"$1"_R2.fastq.gz fastq/"$1"_R2.fastq.gz.og && \
		__ts "Original fastq files are backup up with .og extenion." && \
		__ts "Renaming trimmed files for $1." && \
		mv fastq/"$1"_R1_val_1.fq.gz fastq/"$1"_R1.fastq.gz && \
		mv fastq/"$1"_R2_val_2.fq.gz fastq/"$1"_R2.fastq.gz && \
        mkdir -p fastq/$(dirname $1)/trim_reports && \
        mv fastq/"$1"*_trimming_report.txt fastq/$(dirname $1)/trim_reports

}

# Make function global
export -f __trim_adapt

# Run the qc function
parallel -j $(( $ncpu > 8 ? $(($ncpu/8)) : 1 )) __trim_adapt ::: ${fqlist[@]}

__ts "Adapter trimming completed."

## Align sequences -------------------------------------------------------------

__ts "Aligning reads to reference genome."

# Load modules
module -q reset; module load parallel bwa samtools gatk

# Get fastq samples
fqlist=( $(
	find fastq/ -name "*R?.fastq.gz" | \
		sed 's/^fastq\///; s/_R.\.fastq\.gz$//' | \
		sort -u
) )
__ts "${#fqlist[@]} sequence samples found."

__align_reads(){

    # Get read group strings
    run=( $(zcat fastq/"$1"_R1.fastq.gz | head -n1 | cut -d+ -f1 | tr ":" " ") )
    id=$(echo ${run[2]} ${run[3]} | tr " " ".")
    pu=$(echo ${run[2]} ${run[3]} ${run[9]} | tr " " ".")
    sm=( $(echo $1 | rev | cut -d/ -f1 | rev | tr "_" "-") )


	# Aligning reads to genome
	[[ ! -s align/$1/align.sam && ! -s align/$1/markdup.bam ]] && \
		__ts "Aligning reads to the genome for $1." && \
        mkdir -p align/$1 && \
		bwa mem -t $(( $ncpu > 8 ? 8 : $ncpu )) -v 1 \
            -R "@RG\tID:$id\tPL:ILLUMINA\tPU:$pu\tLB:$sm\tSM:$sm" \
            refs/genome.fna \
            fastq/"$1"_R1.fastq.gz \
            fastq/"$1"_R2.fastq.gz > \
            align/$1/align.sam && \
		__ts "Writing alignment statistics $1." && \
        samtools flagstat align/$1/align.sam > \
            align/$1/alignment_stats.txt || \
		__ts "Alignment file is already present for $1. Skipping alignemnt."
		
	# Sorting and marking the PCR duplicates
	[[ -s align/$1/align.sam && ! -s align/$1/markdup.bam ]] && \
		__ts "Marking duplicates for $1." && \
		gatk MarkDuplicatesSpark -verbosity WARNING --spark-verbosity WARN \
            --remove-sequencing-duplicates \
            -I align/$1/align.sam \
            -O align/$1/markdup.bam \
            -M align/$1/markdup_metrics.txt || \
		__ts "Duplicates are already fixed for $1. Skipping."

	# Removing unsorted BAM file
	[[ -s align/$1/align.sam && -s align/$1/markdup.bam.bai ]] && \
		__ts "Removing unsorted SAM file for $1." && \
		rm align/$1/align.sam		

}

# Make function global
export -f __align_reads

# Run the alignment function
parallel -j $(( $ncpu > 8 ? $(($ncpu/8)) : 1 )) __align_reads ::: ${fqlist[@]}

__ts "Finished alignment."

## Variant calling --------------------------------------------------------

__ts "Variant calling from the alignments."

# Load modules
module -q reset; module load parallel gatk

# Get sorted alignments
alnlist=( $(
	find align/ -name "markdup.bam" | \
		sed 's/^align\///; s/\/markdup\.bam$//'
) )
__ts "${#alnlist[@]} alignments found."

__call_variant(){

	# Variant calling with Haplotypecaller
		[[ ! -s variants/$1/vars.bp.vcf.gz && ! -s variants/$1/vars.bp.vcf ]] && \
		__ts "Variant calling for $1." && \
        mkdir -p variants/$1 && \
        gatk HaplotypeCaller -verbosity WARNING \
            --output-mode EMIT_ALL_ACTIVE_SITES -ERC BP_RESOLUTION \
            -R refs/genome.fna \
            -I align/$1/markdup.bam \
            -O variants/$1/vars.bp.vcf || \
		__ts "Variant file is already present for $1. Skipping."

	# Genotype the variations
	[[ -s variants/$1/vars.bp.vcf && \
        ! -s variants/$1/vars.vcf ]] && \
		__ts "Genotyping $1 variants." && \
        gatk GenotypeGVCFs -verbosity WARNING \
            -R refs/genome.fna \
            -V variants/$1/vars.bp.vcf \
            -O variants/$1/vars.vcf || \
		__ts "$1 is already genotyped. Skipping."

	# Get list of SNPs
	[[ -s variants/$1/vars.vcf && \
        ! -s variants/$1/snps.vcf ]] && \
		__ts "Getting SNPs only for $1." && \
        gatk SelectVariants -verbosity WARNING \
            -select-type SNP \
            -R refs/genome.fna \
            -V variants/$1/vars.vcf \
            -O variants/$1/snps.vcf || \
		__ts "$1 SNP list is already present. Skipping."

}

# Make function global
export -f __call_variant

# Run the polishing function
parallel -j $(( $ncpu > 8 ? $(($ncpu/8)) : 1 )) __call_variant ::: ${alnlist[@]}

__ts "Finished variant calling."

## Getting list of all SNP positions -------------------------------------------

# Load required modules
module -q reset

# Get all SNP only vcfs
snplist=( $(
	find variants/ -name "snps.vcf" | \
    sort
) )
__ts "${#snplist[@]} SNP VCFs found."

__ts "Getting the position of all SNPs."

# Function to index prokaryotic genome
__allsnp_pos() {

    # Get the identifier string for the chromosome
    chr="."
    [[ -s refs/genome.fna ]] && \
        __ts "Getting chromosome ID."
        chr=$(
            awk -vRS=">" -vORS="\n" -vFS="\n" -vOFS="\t" '
                NR>1 {a=$1; $1 = ""; print a, length($0)-NF+1}
            ' refs/genome.fna | \
                sort -k2n | \
                head -n1 | \
                cut -f1 -d" "
        ) || \
        __ts "Unable to generate chromosome ID from reference genome."

    # Generating filter table
    [[ ! -s variants/filter.table ]] && \
        __ts "Generating SNP filter list." && \
        cat ("$@") | \
        grep -v "^#" | \
        grep $chr | \
        awk -vqd=$qd -vadp=$adp '{OFS="\t";
            split($8, b, ";");
            for (j in b) {split(b[j], c, "="); d[c[1]] = c[2]};
            split($10,e,":"); split(e[2],f,",");
            if (d["QD"] >= qd && f[2] >= adp && f[1] < adp) print $1, $2;
        }' | \
        sort -k1,1 -k2,2n | uniq > \
        variants/filter.table || \
        __ts "SNP filter list is already present. Skipping."

}

# Running the SNP allele table generation
__allsnp_pos $snplist

## Generate allele table -------------------------------------------------------

# Load required modules
module -q reset

__ts "Generating allele tables."

# Get list of BP variants
varlist=( $(
    find variants/ -name "vars.bp.vcf" | \
		sed 's/^variants\///; s/\/vars\.bp\.vcf$//' | \
        sort -u
) )
__ts "${#varlist[@]} variants found."

# Function to generate allele table
__allele_table() {

    # Preparing allele table
    [[ ! -s ariants/$1/allele.table ]] && \
        __ts "Preparing allele table for $1." && \
        awk -vOFS="\t" -vdp=$dp -vrdp=$rdp '
            BEGIN{print "CHROM", "POS", "REF", "ALT", "QD", "RDP", "ADP"}
            NR==FNR{a[$1 FS $2]=$0; next}
            ($1 FS $2) in a{
                split($5, b, ",")
                split($10, c, ":")
                split(c[2], d, ",")
                if (length(c)<5 || length($4)>1) {print $1, $2, ".", ".", 0, 0, 0; next}
                for (i=1; i<=length(b); i++) e[b[i]]=d[i+1]
                alt=b[1]
                for (f in b) alt=length(alt)<=length(b[f])?alt:b[f]
                if (length(alt)>1 && alt!="<NON_REF>") {print $1, $2, $4, ".", 0, 0, 0; next}
                print $1, $2, $4, c[3]<=dp?".":alt=="<NON_REF>"?$4:c[3]-e[alt]>=rdp?$4:alt, c[3]<=dp?"0":$6/c[3], d[1], e[alt]
            }
        ' variants/filter.table variants/$1/vars.bp.vcf > \
            variants/$1/allele.table || \
        __ts "Allele table is already present for $1. Skipping."

}

# Make function global
export -f __allele_table

# Run the allele tabulation function
parallel -j $(( $ncpu > 8 ? 8 : $ncpu )) __allele_table ::: ${varlist[@]}

__ts "Completed tabulating alleles."

## Generate SNP only multiFASTA ------------------------------------------------

# Load required modules
module -q reset

# Getting list of allele tables
tablist=( $(
    find variants/ -name "allele.table" | \
        sort
) )
__ts "${#tablist[@]} allele tables found."

# Function to generate combined FASTA
__combn_snps() {

    # Getting list of tables from args
    list=("$@")

    # Generating combined allele table
    [[ ! -s variants/combined.table ]] && \
        __ts "Generating combined allele table."
        echo -e "CHROM\tPOS\t"${list[*]} | \
            tr " " "\t" | \
            sed 's/variants\///g; s/\/allele\\.table//g' > \
            variants/combined.table && \
        paste $(printf ' %q' "${list[@]}") -d"\t" | \
            awk -vadp=$adp -vqd=$qd '
                NR>1{
                    printf "%s\t%s\t", $1, $2;
                    for(i=0;i<=NF;i+=7) printf "%s\t%s", ($(i+5)>=qd && $(i+7)>adp && $(i+6)<adp)?$(i+4):$(i+3), (i+7>NF?"\n":FS)
                }
            ' >> variants/combined.table || \
        __ts "Unable to generate chromosome ID from reference genome."

    # Generating FASTA from allele table
    [[ ! -s variants/combined.fasta ]] && \
        __ts "Generating combined FASTA." && \
        cut -f1,2 --complement variants/combined.table | \
            grep -v "\." | \
            grep -v "\*" | \
            awk '{for(c=1;c<=NF;c++){a[c,NR]=$c}if(max_nf<NF){max_nf=NF}}END{for(r=1;r<=max_nf;r++){for(c=1;c<=NR;c++){printf("%s ", a[r, c])}print ""}}' | \
            awk '{split($1, a, "_"); $1=""; gsub(" ", ""); print ">"a[1]"\n"$0}' > \
            variants/combined.fasta || \
        __ts "Combined fasta is already present. Skipping."

}

# Running the FASTA generation function
__combn_snps ${tablist[*]}

## Generate phylogenetic tree --------------------------------------------------

# Load required modules
module -q reset; module load raxml newick_utils/1.6

__ts "Generating phylogenetic tree."

# Function to generate tree
__tree_raxml() {

    # Generating phytogenetic tree
    [[ ! -s output/RAxML_bootstrap.combined.tree && \
        ! -s output/RAxML_result.combined.tree ]] && \
        __ts "Generating the tree."
        raxmlHPC -d -p 144 -b 144 -N $boot -m GTRCATI \
            -s variants/combined.fasta \
            -n combined.tree \
            -w output || \
        __ts "Phylogenetic tree is already present. Skipping."

    # Display the tree (if bootstraped)
    [[ -s output/RAxML_bootstrap.combined.tree ]] && \
        nw_display output/RAxML_bootstrap.combined.chr.tree2

    # Display the tree
    [[ -s output/RAxML_result.combined.tree ]] && \
        nw_display output/RAxML_result.combined.chr.tree2

}

# Running the phylogenetic tree function
__tree_raxml
