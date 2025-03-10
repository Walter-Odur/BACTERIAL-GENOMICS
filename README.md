<h1 align="center">
    <b>VARIANT CALLING PIPELINE</b> <br>
    <b>DOCUMENTATION</b>
</h1>

<h3 align="center">
    Author: [Your Name] <br>
    Affiliation: [Your Institution]
</h3>

## **Background**
<div align="justify">
Variant calling is the process of identifying genetic variations such as 
single nucleotide polymorphisms (SNPs) and insertions/deletions (INDELs) from 
next-generation sequencing (NGS) data. This pipeline provides a streamlined 
workflow for processing raw sequencing data, aligning reads to a reference genome, 
and identifying variants using bioinformatics tools such as BWA, Samtools, and BCFtools.  
</div>

## **Pipeline Overview**
The pipeline consists of the following major steps:
1. **Quality Control**: Assessing raw read quality using FastQC.
2. **Read Trimming**: Removing low-quality bases and adapters using Fastp.
3. **Alignment**: Mapping reads to a reference genome using BWA-MEM.
4. **BAM Processing**: Sorting, indexing, and converting SAM to BAM using Samtools.
5. **Variant Calling**: Generating a VCF file using BCFtools.

---

## **📌 Prerequisites**
Ensure you have the following tools installed before running the pipeline:

- [BWA](http://bio-bwa.sourceforge.net/) - Burrows-Wheeler Aligner for read mapping.
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - Quality control of sequencing data.
- [Fastp](https://github.com/OpenGene/fastp) - Trimming of low-quality reads.
- [Samtools](http://www.htslib.org/) - Processing of BAM/SAM files.
- [BCFtools](http://www.htslib.org/) - Variant calling and filtering.

💡 **Note:** This pipeline assumes you are working with **paired-end sequencing data** stored in a `reads/` directory.

---

## **The Pipeline Workflow**
<figure>
<img src="variant_calling_pipeline.png" width="732" 
alt="This figure represents the step-by-step workflow of the variant calling pipeline." />
<figcaption aria-hidden="true">Workflow from raw reads to variant identification.</figcaption>
</figure>

# **The Bash Pipeline: Variant Calling from Raw Reads**

- The **Bash script** (`variant_calling.sh`) for running the complete pipeline is provided
  [**HERE**](https://github.com/your-repo).  
  Below is a breakdown of all steps included in the script.

---

```bash
# Index the reference genome
bwa index ref/GCF_000195955.2_ASM19595v2_genomic.fna

# Quality control with FastQC
fastqc -o Results/fastQC reads/ERR987696_1.fastq reads/ERR987696_2.fastq

# Read trimming with Fastp
fastp -i reads/ERR987696_1.fastq -I reads/ERR987696_2.fastq \
      -o Results/trim/ERR987696_1_trim.fastq -O Results/trim/ERR987696_2_trim.fastq

# Read alignment with BWA
bwa mem -t 16 -aM ref/GCF_000195955.2_ASM19595v2_genomic.fna \
        Results/trim/ERR987696_1_trim.fastq Results/trim/ERR987696_2_trim.fastq \
        > Results/bams/ERR987696.sam

# Convert SAM to BAM
samtools view -bS Results/bams/ERR987696.sam > Results/bams/ERR987696.bam

# Sort and index BAM files
samtools sort Results/bams/ERR987696.bam -o Results/bams/ERR987696_sorted.bam
samtools index Results/bams/ERR987696_sorted.bam

# Generate mpileup file
bcftools mpileup -f ref/GCF_000195955.2_ASM19595v2_genomic.fna \
                 Results/bams/ERR987696_sorted.bam > Results/mpileups/ERR987696.mpileup

# Perform variant calling
bcftools call -v -m -O v -o Results/vcfs/ERR987696_variants.vcf Results/mpileups/ERR987696.mpileup
