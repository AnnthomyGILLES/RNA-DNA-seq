#!/bin/bash

"
This is a standard GATK and lofreq variant call pipeline with snpEff annotations.

We begin by mapping the sequence reads to the reference genome to produce a file in SAM/BAM format sorted by coordinate. Next, we mark duplicates to mitigate biases introduced by data generation steps such as PCR amplification. Then we perform local realignment around indels, because the algorithms that are used in the initial mapping step tend to produce various types of artifacts in the regions around indels. Finally, we recalibrate the base quality scores, because the variant calling algorithms rely heavily on the quality scores assigned to the individual base calls in each sequence read.

Raw data needs to be labeled with the exact number of characters as in the example bellow:

SN7640335_40788_E9_sequence.fq
" > /dev/null 2>&1

adapters_file=/beegfs/common/genomes/adapters/TruSeqAdapters.txt # TruSeq Adapters
oriFasta=/beegfs/group_bit_gu/data/projects/departments/Bioinformatics/Ann_pipeline/Variant_calling_Dnav2/scripts/Caenorhabditis_elegans.WS220.61.dna.toplevel.fa
faref=reference.fa # name to be attributed to the reference fasta on your project folder
fadic=reference.dict # name to be attributed to the reference dictionary on your project folder
pic=/beegfs/common/software/picard-tools/1.137/picard.jar
#vcf=/beegfs/group_bit_gu/data/projects/departments/Bioinformatics/Ann_pipeline/Variant_calling_Dna/Haplotype_Caller/output_filtred.vcf 
vcf=
snpEffdb=/beegfs/group_bit_gu/home/AGilles/Ann_pipeline/Variant_calling_YiDR/snpEff/WBcel235_82/

# Create output folders
mkdir ../fastqc_output > /dev/null 2>&1
mkdir ../raw_data > /dev/null 2>&1
mkdir ../flexbar_output > /dev/null 2>&1
mkdir ../slurm_logs > /dev/null 2>&1
mkdir ../tmp > /dev/null 2>&1
mkdir ../replace_groups > /dev/null 2>&1
mkdir ../sortsam_output > /dev/null 2>&1
mkdir ../markduplicates_output > /dev/null 2>&1
mkdir ../snpeff_output > /dev/null 2>&1
mkdir ../others > /dev/null 2>&1
mkdir ../realign_output > /dev/null 2>&1
mkdir ../recalibrate_output > /dev/null 2>&1
mkdir ../variantscall_output > /dev/null 2>&1
mkdir ../variantscall_output/gvcfs> /dev/null 2>&1
mkdir ../variant_filtration > /dev/null 2>&1
mkdir ../variant_filtration/selectvariants > /dev/null 2>&1
mkdir ../variant_filtration/catvariants > /dev/null 2>&1
mkdir ../variant_filtration/variantfiltration > /dev/null 2>&1
mkdir ../lofreq_output > /dev/null 2>&1
mkdir ../Haplotype_Caller_output > /dev/null 2>&1
mkdir ../join_genotype > /dev/null 2>&1
mkdir ../bwa_output > /dev/null 2>&1
mkdir ../others/f_alt_output > /dev/null 2>&1


# Paths
top=$(readlink -f ../)/
raw=/beegfs/group_bit_gu/home/AGilles/Ann_pipeline/Variant_calling_Dna/raw_data/
filt=$(readlink -f ../variant_filtration)/
qc=$(readlink -f ../fastqc_output)/
snpe=$(readlink -f ../snpeff_output)/
logs=$(readlink -f ../slurm_logs)/
tmp=$(readlink -f ../tmp)/
fb=$(readlink -f ../flexbar_output)/
bm=$(readlink -f ../bwa_output)/
rg=$(readlink -f ../replace_groups)/
ss=$(readlink -f ../sortsam_output)/
md=$(readlink -f ../markduplicates_output)/
others=$(readlink -f ../others)/
real=$(readlink -f ../realign_output)/
rec=$(readlink -f ../recalibrate_output)/
var=$(readlink -f ../variantscall_output)/
lofreq=$(readlink -f ../lofreq_output)/
HC=$(readlink -f ../Haplotype_Caller_output)/
jg=$(readlink -f ../join_genotype)/

f_alt=$(readlink -f ../others/f_alt_output)/
gvcfs=$(readlink -f ../variantscall_output/gvcfs)/
selVar=$(readlink -f ../variant_filtration/selectvariants)/
varFilt=$(readlink -f ../variant_filtration/variantfiltration)/
catVar=$(readlink -f ../variant_filtration/catvariants)/
# Load required software
module load vcftools
module load IGV
module load Java
module load SAMtools
module load STAR
module load FastQC
module load Flexbar
module load R
module load lofreq
module load BWA 
module load tabix
module load BCFtools

snpEff=/beegfs/group_bit_gu/data/projects/departments/Bioinformatics/Ann_pipeline/Variant_calling_YiDR/snpEff/snpEff.jar
GA=/beegfs/common/software/GATK/3.4-46/GenomeAnalysisTK.jar

# Generate indexes
printf '\nGenerating indexes\n\n'  
mkdir ../others > /dev/null 2&>1
cd ../others

srun cp ${oriFasta} ${faref}
srun --cpus-per-task=18 bwa index ${faref} 
srun java -jar ${pic} CreateSequenceDictionary R=${faref} O=${fadic}
srun samtools faidx ${faref}



fas=$(ls *.*)
faref=$(readlink -f ${faref})
IDS=
id1=
id=
cd ${raw}

for f in $(ls *_sequence.fq);do 
    echo "#!/bin/bash
    
    printf '\nFastQC\n\n'
    fastqc -t 4 -o ${qc} ${f}
    
    printf '\nFlexbar\n\n'
    flexbar -r ${f} -t ${fb}${f::(-3)} -n 18 -a ${adapters_file} -ao 10 -u 5 -q 20 -f i1.8 -ae ANY

    printf '\nMap to reference\n\n'
    #bwa mem -M -t 18 ${ref} ${fb}${f::(-3)} > ${bm}${f::(-12)}.sam  
    bwa mem -M -t 18 -R '@RG\tID:${f:16:(-12)}\tSM:${f:16:(-12)}\tPL:illumina\tLB:lib1\tPU:unit1' ${faref} ${fb}${f::(-3)}.fastq > ${bm}${f::(-12)}.sam 
 
    printf '\nSort SAM\n\n'
	java -jar ${pic} SortSam SORT_ORDER=coordinate INPUT=${bm}${f::(-12)}.sam OUTPUT=${ss}${f::(-12)}.bam TMP_DIR=`pwd`/tmp

	printf '\nRead groups\n\n' 
    java -jar ${pic} AddOrReplaceReadGroups I=${ss}${f::(-12)}.bam O=${rg}${f::(-12)}.bam SO=coordinate RGID=${f:16:(-12)} RGLB=lib1 RGPL=illumina RGPU=unit1 \
RGSM=${f:16:(-12)}  

    printf '\nMark duplicates\n\n' 
    java -jar ${pic} MarkDuplicates I=${rg}${f::(-12)}.bam O=${md}${f::(-12)}.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${md}metrics.txt
    cd ${md}
    
    printf '\nIndexing bam\n\n' 
    java -jar ${pic} BuildBamIndex INPUT=${f::(-12)}.bam

	if [ -z ${vcf} ];then
    	printf '\nHaplotype Caller to generate G.vcf\n\n'  
		java -jar ${GA} -R ${faref} -T HaplotypeCaller -I ${md}${f::(-12)}.bam -ERC GVCF -stand_call_conf 20 -stand_emit_conf 20 -o ${HC}${f::(-12)}.indels.g.vcf
    fi
    rm ${tmp}${f}_p1.sh
" > ${tmp}${f}_p1.sh
    chmod 755 ${tmp}${f}_p1.sh
    rm ${logs}${f}_p1.*.out > /dev/null 2>&1
    ID1=$(sbatch -o ${logs}${f}_p1.%j.out --cpus-per-task=18 ${tmp}${f}_p1.sh)
	id1=${id1}:${ID1:20}  
	echo ${ID1}
done 
	echo "#!/bin/bash
	cd ${HC}
	if [ -z ${vcf} ];then
		ls *.g.vcf > mygvcfs.list

		printf '\nGenotypeGVCFs to generate G.vcf\n\n'  
		java -jar ${GA} -T GenotypeGVCFs -R ${faref} --variant mygvcfs.list -o results_GenotypeGVCFs.vcf

		printf '\VariantFiltration to generate G.vcf\n\n'  
		java -jar ${GA} -T VariantFiltration -R ${faref} -V results_GenotypeGVCFs.vcf -window 35 -cluster 3 -filterName FS -filter 'FS > 30.0' -filterName QD \
	-filter 'QD < 2.0' -o output_filtred.vcf 

	fi

	# Generate indexes
	printf '\nGenerating indexes\n\n'  
	cd ${others}

	srun vcftools --vcf ${HC}output_filtred.vcf  --out output_filtred_indels --recode --recode-INFO-all --keep-only-indels


	srun igvtools index output_filtred_indels.recode.vcf

	
	rm ${tmp}Merge_p2.sh
    " > ${tmp}Merge_p2.sh
    chmod 755 ${tmp}Merge_p2.sh
    rm ${logs}Merge_p2.*.out > /dev/null 2>&1
    id=$(sbatch -d afterok${id1} -o ${logs}Merge_p2.%j.out ${tmp}Merge_p2.sh)
    id2=${id2}:${id:20}
    echo ${id}

cd ${raw}
for f in $(ls  *_sequence.fq);do
    echo "#!/bin/bash
	deletions=$(readlink -f output_filtred_indels.recode.vcf)
	vcf=$(readlink -f ../Haplotype_Caller_output/output_filtred.vcf)
    printf '\nRealigner Target Creator\n\n' 
	java -jar ${GA} -nt 18 -T RealignerTargetCreator -known ${deletions} -R ${faref} -I ${md}${f::(-12)}.bam -o ${others}${f::(-12)}.list 

    # IndelRealigner has a marginal effect on the number of indels which are recovered,
    # it is memory greedy and therefore it tends to fail and it can be skipped.  
    printf '\nIndel Realigner\n\n' 
    java -jar ${GA} -T IndelRealigner -known ${deletions} -R ${faref} -I ${md}${f::(-12)}.bam -targetIntervals ${others}${f::(-12)}.list \
-o ${real}${f::(-12)}.bam
    ####rm ${real}${f::(-12)}*
    ####cp ${md}${f::(-12)}.bam ${real}${f::(-12)}.bam
    cd ${real}
    printf '\nBuild Bam index\n\n'
    java -jar ${pic} BuildBamIndex INPUT=${f::(-12)}.bam

    printf '\nFirst Base Recalibrator\n\n'
    java -Xmx60g -jar ${GA} -T BaseRecalibrator -lowMemoryMode -R ${faref} -I ${real}${f::(-12)}.bam -knownSites ${vcf} -o ${others}${f::(-12)}.table

    printf '\nSecond Base Recalibrator\n\n'     
    java -jar ${GA} -T BaseRecalibrator -R ${faref} -I ${real}${f::(-12)}.bam -knownSites ${vcf} -BQSR ${others}${f::(-12)}.table \
-o ${others}${f::(-12)}_post.table

    printf '\nAnalyze Covariates\n\n'    
    java -jar ${GA} -T AnalyzeCovariates -R ${faref}  -l DEBUG -before ${others}${f::(-12)}.table -after ${others}${f::(-12)}_post.table \
-plots ${others}${f::(-12)}.recal_plots.pdf

    printf '\nPrint Reads\n\n'  
    java -jar ${GA} -T PrintReads -R ${faref} -I ${real}${f::(-12)}.bam -BQSR ${others}${f::(-12)}.table -o ${rec}${f::(-12)}.bam 

    printf '\nHaplotype Caller\n\n'  
	java -jar ${GA} -T HaplotypeCaller -R ${faref} -I ${rec}${f::(-12)}.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 \
-ERC GVCF -o ${gvcfs}${f::(-12)}.g.vcf -bamout ${var}${f::(-12)}.bam --variant_index_type LINEAR --variant_index_parameter 128000

	printf '\nGenotypeGVCFs\n\n'  
	java -jar ${GA} -T GenotypeGVCFs -R ${faref} --variant ${gvcfs}${f::(-12)}.g.vcf -o  ${jg}${f::(-12)}.vcf


	###############			HARD FILTERS			#####################
	cd ${others}
    printf '\nHard Filters: SelectVariants\n\n'  
 	java -jar ${GA} -T SelectVariants -R ${faref} -V ${jg}${f::(-12)}.vcf -selectType SNP -o ${var}${f:16:(-12)}_raw_snps.vcf
    java -jar ${GA} -T SelectVariants -R ${faref} -V ${jg}${f::(-12)}.vcf -selectType INDEL -o ${var}${f:16:(-12)}_raw_INDELS.vcf   
    
	printf '\nHard Filters: VariantFiltration\n\n'  
	#java -jar ${GA} -T VariantFiltration -R ${faref} -V ${var}${f:16:(-12)}_raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName '${f:16:(-12)}_snps' -o ${varFilt}${f:16:(-12)}_filt_snps.vcf  
	#java -jar ${GA} -T VariantFiltration -R ${faref} -V ${var}${f:16:(-12)}_raw_INDELS.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --filterName '${f:16:(-12)}_indels' -o ${varFilt}${f:16:(-12)}_filt_INDELS.vcf

	printf '\nHard Filters: SelectVariants\n\n'  
	java -jar ${GA} -T SelectVariants -R ${faref} --variant ${var}${f:16:(-12)}_raw_snps.vcf -select 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' -o ${selVar}${f:16:(-12)}_PASS_snps.vcf
	java -jar ${GA} -T SelectVariants -R ${faref} --variant ${var}${f:16:(-12)}_raw_INDELS.vcf -select 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' -o ${selVar}${f:16:(-12)}_PASS_INDELS.vcf

	java -cp ${GA} org.broadinstitute.gatk.tools.CatVariants -R ${faref} -V ${selVar}${f:16:(-12)}_PASS_snps.vcf -V ${selVar}${f:16:(-12)}_PASS_INDELS.vcf -out ${catVar}${f:16:(-12)}.pass.vcf -assumeSorted


	###############			SNPEFF ANNOTATION			#####################

	java -Xmx4g -jar ${snpEff} -ud 0 -v ws220.64 ${catVar}${f:16:(-12)}.pass.vcf > ${snpe}${f:16:(-12)}.pass.ann.vcf

	###############			LOFRED ANNOTATION			#####################

    printf '\nLoFreq\n'
    lofreq call --call-indels --no-default-filter --verbose -f ${faref} -o ${lofreq}${f::(-12)}.vcf ${rec}${f::(-12)}.bam

    printf '\nsnpEff on LoFreq output\n'
    java -Xmx4g -jar ${snpEff} -ud 0 -v ws220.64 ${lofreq}${f::(-12)}.vcf > ${lofreq}${f:16:(-12)}.ann.vcf

	java -jar ${GA} -T FastaAlternateReferenceMaker -R ${faref} -o ${f_alt}${f:16:(-12)}.g.fa -V ${catVar}${f:16:(-12)}.pass.vcf

    #printf '\VariantFiltration\n\n'  
	#java -jar ${GA} -T VariantFiltration -R ${faref} -V ${jg}${f::(-12)}.vcf -window 35 -cluster 3 -filterName FS -filter 'FS > 30.0' -filterName QD \
-filter 'QD < 2.0' -o ${filt}${f::(-12)}.vcf 
    
    rm ${tmp}${f}_p3.sh
    " > ${tmp}${f}_p3.sh
    chmod 755 ${tmp}${f}_p3.sh
    rm ${logs}${f}_p3.*.out > /dev/null 2>&1
    id3=$(sbatch -d afterok${id2} -o ${logs}${f}_p3.%j.out ${tmp}${f}_p3.sh)
	echo ${id3}
done

exit

