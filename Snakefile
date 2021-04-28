#################################
#snakefile for implementing PacBio variant calling pipeline
#Justin Lack, justin.lack@nih.gov
#NIAID Collaborative Bioinformatics Resource (NCBR)
#July 19, 2019
################################
################################

DATA_SAMPLES = ["Thvx04","Thvx32"]
#data_path="/data/NCBR/rawdata/NCBR-118/round4_Jan10_2021/"
refs_path="/data/NCBR/rawdata/NCBR-118/"
ref_fasta = refs_path+"Homo_sapiens_assembly38.fasta"
gff_file = refs_path+"transcriptome.gff3.gz"
configfile:"references.json"

rule all:
    input:
        expand("{samples}/{samples}.R1.bam", samples=DATA_SAMPLES),
        expand("{samples}/{samples}.R2.bam", samples=DATA_SAMPLES),
        expand("{samples}/{samples}.insertions.bed", samples=DATA_SAMPLES),
        expand("{samples}/QC/{samples}HIV_further_filtered_R2.bamtools", samples=DATA_SAMPLES),
        "merged_sorted_insertions.annotations.txt"
    params: rname="final"
    output: "multiqc_report.html"
    shell: """module load multiqc/1.8; multiqc ."""

rule merge:
    input: dir="rawdata/{samples}",
    output: r1bam = "{samples}/{samples}.R1.bam",
            r2bam = "{samples}/{samples}.R2.bam",
    threads:8
    params: sample = "{samples}", rname="merge"
    shell: """
           mkdir -p {params.sample}
           ls rawdata/{params.sample}/*R1*.molbar.trimmed.deduped.bam > {params.sample}/r1_files
           ls rawdata/{params.sample}/*R2*.molbar.trimmed.deduped.bam > {params.sample}/r2_files
           module load samtools
           samtools merge {output.r1bam} -b {params.sample}/r1_files
           samtools merge {output.r2bam} -b {params.sample}/r2_files
           """

rule process:
    input: r1bam = "{samples}/{samples}.R1.bam",
           r2bam = "{samples}/{samples}.R2.bam",
    output: bed = "{samples}/{samples}.insertions.bed",
    params: sample = "{samples}", rname="process"
    shell: """
           module load python/2.7
           module load samtools
           python NCBR_hiv_integration/pairwise2_v0.5.3_RC_compatible.py {input.r2bam} -v
           samtools view -h -o {params.sample}/{params.sample}.R1.sam {input.r1bam}
           python NCBR_hiv_integration/select_pairs_with_GT24bp_hg19_mapping.py {params.sample}/{params.sample}.R1.sam {params.sample}/{params.sample}HIV_aligned.sam
           mv {params.sample}HIV_further_filtered_R2.sam {params.sample}/{params.sample}HIV_further_filtered_R2.sam
           samtools view -bS {params.sample}/{params.sample}HIV_further_filtered_R2.sam > {params.sample}/{params.sample}HIV_further_filtered_R2.bam
           samtools sort {params.sample}/{params.sample}HIV_further_filtered_R2.bam > {params.sample}/{params.sample}HIV_further_filtered_R2.sorted.bam
           samtools index {params.sample}/{params.sample}HIV_further_filtered_R2.sorted.bam
           module load bedtools/2.29.2
           bamToBed -i {params.sample}/{params.sample}HIV_further_filtered_R2.sorted.bam > {params.sample}/{params.sample}.insertions.bed
           """

rule qc:
    input: r1bam = "{samples}/{samples}.R1.bam",
           r2bam = "{samples}/{samples}.R2.bam",
           bed = "{samples}/{samples}.insertions.bed",
    output: stats="{samples}/QC/{samples}HIV_further_filtered_R2.bamtools"
    params: sample = "{samples}", adapters=config['references']['fastqc_adapters'], rname="qc"
    shell: """
           mkdir -p {params.sample}/QC
           module load fastqc/0.11.8
           fastqc -o {params.sample}/QC --threads 2 -f bam --contaminants {params.adapters} {input.r1bam}
           fastqc -o {params.sample}/QC --threads 2 -f bam --contaminants {params.adapters} {input.r2bam}
           module load samtools/1.9
           samtools flagstat {input.r1bam} > {params.sample}/QC/{params.sample}.R1.flagstat
           samtools flagstat {input.r2bam} > {params.sample}/QC/{params.sample}.R2.flagstat
           samtools flagstat {params.sample}/{params.sample}HIV_further_filtered_R2.bam > {params.sample}/QC/{params.sample}HIV_further_filtered_R2.flagstat
           module load bamtools/2.5.1
           bamtools stats -in {input.r1bam} > {params.sample}/QC/{params.sample}.R1.bamtools
           bamtools stats -in {input.r2bam} > {params.sample}/QC/{params.sample}.R2.bamtools
           bamtools stats -in {params.sample}/{params.sample}HIV_further_filtered_R2.bam > {params.sample}/QC/{params.sample}HIV_further_filtered_R2.bamtools
           """

rule prep_bed:
    input: bed = "{samples}/{samples}.insertions.bed",
    output: bed="{samples}/{samples}.insertions.collapsed.bed",
    params: sample = "{samples}",rname="prep"
    shell: """
           module load perl
           perl NCBR_hiv_integration/preprocess_insertions2.pl {input.bed} {output.bed} {params.sample}
           """

rule annotation:
    input: expand("{samples}/{samples}.insertions.collapsed.bed", samples=DATA_SAMPLES),
    output: vep="merged_sorted_insertions.annotations.txt",
    params: rname="vep"
    shell: """
           cat {input} > merged_insertions.bed
           module load bedtools/2.29.2
           bedtools sort -i merged_insertions.bed > merged_sorted_insertions.bed
           module load perl
           perl /data/NCBR/rawdata/NCBR-118/scripts/alt_make_vep_input.pl merged_sorted_insertions.bed merged_sorted_insertions_vep.bed
           module load VEP/100
           vep -i merged_sorted_insertions_vep.bed -o {output.vep} --species human --everything --assembly GRCh37 --offline --cache --dir_cache $VEP_CACHEDIR --tab --fasta $VEP_CACHEDIR/GRCh37.fa --pick
           """
