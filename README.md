<h1 align="center">
    <b>BACTERIAL GENOMICS</b> <br>
    <b>VARIANT CALLING PIPELINE DOCUMENTATION</b>
</h1>

<h3 align="center">
    Author: Walter Odur <br>
    Affiliation: ACE-Uganda
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
- [BCFtools](https://samtools.github.io/bcftools/bcftools.html) - Variant calling and filtering.

💡 **Note:** This pipeline assumes you are working with **paired-end sequencing data** stored in a `reads/` directory.

---

## **The Pipeline Workflow**
<figure>
<img src="variant_calling_pipeline.png" width="732" 
alt="This figure represents the step-by-step workflow of the variant calling pipeline." />
<figcaption aria-hidden="true">Workflow from raw reads to variant identification.</figcaption>
</figure>

# **The Detailed Steps**

- The **Bash script** (`variant_calling.sh`) for running the complete pipeline is provided
  [**HERE**](https://github.com/Walter-Odur/BACTERIAL-GENOMICS).  
  Below is a breakdown of all steps included in the script.

---

## **Step 1: Index the Reference Genome**
<div align="justify">
Indexing the reference genome is essential for efficient read alignment.
BWA uses indexed files to rapidly find alignment positions without scanning
the entire genome.
</div>

```bash
bwa index ref/GCF_000195955.2_ASM19595v2_genomic.fna
