#!/bin/bash


# Submit the slurm script

#SBATCH --job-name=VariantCalling
#SBATCH --output=var_output.txt
#SBATCH --error=var_errors.txt
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=16


outdir=Results
mkdir -p ${outdir}/{fastQC,trim,bams,vcfs,mpileups}

ls -1 reads/*.fastq |cut -d '_' -f1 |sort -u |xargs -i{} basename {} >samples.txt

samples=$(cat samples.txt)

# Reference file was obtained using the following link:
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz

ref=$(ls ref/*.{fa,fasta,fna} 2>/dev/null)


# STEP1: INDEX THE REFERENCE
if ! ls  ref/*.bwt 1>/dev/null 2>&1
then
                echo "Indexing the reference"
                bwa index $ref
else
                echo "reference already indexed"
fi


# STEP2: QUALITY CONTROL

for sample in ${samples[@]}
do
        if ! ls ${outdir}/fastQC/${sample}*_fastqc.zip 1>/dev/null 2>&1
        then
                        echo "Running fastqc on ${sample}"
                        fastqc \
                                -o ${outdir}/fastQC \
                                reads/${sample}*.fastq

        else
                        echo "fastqc already done for ${sample}"
        fi
done

# STEP3: MULTIQC
# multiqc \
#         -o ${outdir}/fastQC \
#         ${outdir}/fastQC


# STEP3: TRIMMING POOR QUALITY READS

for sample in ${samples[@]}
do
        if ! ls ${outdir}/trim/${sample}*trim.fastq 1>/dev/null 2>&1
        then
                        echo "trimming sample: ${sample}"
                        fastp \
                                -i reads/${sample}_*1.fastq \
                                -I reads/${sample}_*2.fastq \
                                -o ${outdir}/trim/${sample}_1_trim.fastq \
                                -O ${outdir}/trim/${sample}_2_trim.fastq
        else
                        echo "${sample} already trimmed"
        fi
done

# STEP4: ALIGNING READS TO THE REFERENCE GENOME

for sample in ${samples[@]}
do
        if ! ls ${outdir}/bams/${sample}*.bam 1>/dev/null 2>&1
        then

                        echo "Aligning ${sample} to the reference genome"
                        bwa mem \
                                -t 16 -aM $ref -R "@RG\tID:$sample\tSM:$sample" \
                                ${outdir}/trim/${sample}_1_trim.fastq \
                                ${outdir}/trim/${sample}_2_trim.fastq | \
                                # convert sam to bam
                                samtools view -bS - > ${outdir}/bams/${sample}.bam
        else
                        echo "${sample} already aligned"

        fi
done


# STEP5: SORT AND INDEX THE BAM

for sample in ${samples[@]}
do
        if ! ls ${outdir}/bams/${sample}*_sorted.bam.bai 1>/dev/null 2>&1
        then
                        echo "sorting and indexing ${sample}.bam"
                        samtools sort ${outdir}/bams/${sample}.bam \
                                > ${outdir}/bams/${sample}_sorted.bam

                        samtools index \
                                ${outdir}/bams/${sample}_sorted.bam
        else
                        echo "${sample}.bam already sorted and indexed"
        fi
done

# Create a list of bam files
ls ${outdir}/bams/*sorted.bam > bams.txt

variantSampleCounts=$(grep '^#CHROM' "${outdir}/vcfs/variants.vcf" | cut -f10- | tr '\t' '\n' | wc -l)
pileupSampleCounts=$(grep '^#CHROM' "${outdir}/mpileups/mpileup.vcf" | cut -f10- | tr '\t' '\n' | wc -l)
sampleCounts=$(wc -l < samples.txt)

# STEP6: GENERATE MPILEUP FILES

if [ "${pileupSampleCounts}" -ne "${sampleCounts}" ]
then
                echo "generating mpileup files for all samples"
                bcftools mpileup \
                        -o ${outdir}/mpileups/mpileup.vcf \
                        -f ${ref} \
                        --bam-list bams.txt

else
                echo "mpileup files already generated"
fi


# STEP7: CALLING VARIANTS

if [ ${variantSampleCounts} -ne ${sampleCounts} ]
then
                echo "Performing variant calling on all samples"
                bcftools call \
                        -v -m -O v \
                        -o ${outdir}/vcfs/variants.vcf \
                        ${outdir}/mpileups/mpileup.vcf
else
                echo "Variants already called"
fi

echo "Everything successfully completed"
