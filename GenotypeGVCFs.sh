#! /bin/sh

# Combine g.VCF files 
# Institute: Mayo Clinic


# reference genome (b37), dbsnp
gatk_db=/data5/dom/rs-cardio/kullo/emerge2/dataBase/public/resource_bundle
refGenome=$gatk_db/human_g1k_v37_decoy.fasta
dbsnp=$gatk_db/dbsnp_138.b37.vcf
hapmap=$gatk_db/hapmap_3.3.b37.vcf
omni=$gatk_db/1000G_omni2.5.b37.vcf
G_1K=$gatk_db/1000G_phase1.indels.b37.vcf
mills=$gatk_db/Mills_and_1000G_gold_standard.indels.b37.vcf

# programs: bwa (ver: 0.7.13), GATK (ver: 3.5), and picard (2.8.0)
software_dir=/data5/dom/rs-cardio/kullo/emerge2/softWare
GATK=$software_dir/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar
#GATK=$software_dir/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
OUT_DIR=/data5/dom/rs-cardio/kullo/emerge2/users/m164988/PADseq_mine



# # # # # # # # # # # # #
# # # # functions # # # #
# # # # # # # # # # # # #

combineGVCFs () {
	chr=$1
	vcf_list=`find -name "*_${chr}.g.vcf" -exec echo "--variant " {} \;`
	java -Xmx128G -jar $GATK -T CombineGVCFs -R $refGenome $vcf_list --disable_auto_index_creation_and_locking_when_reading_rods -o cohort_chr${chr}.g.vcf
}

genotypeGVCFs () {
	chr=$1
	java -Xmx128G -jar $GATK -T GenotypeGVCFs -R $refGenome --dbsnp $dbsnp --variant cohort_chr${chr}.g.vcf -o chr${chr}.vcf 
}

selectVariants () {
	d=$1
	# select variants according to type	
	java -Xmx96G -jar $GATK \
		-T SelectVariants \
		-R $refGenome \
		-V chr${d}.vcf \
		-o chr${d}_snvs.vcf \
		-selectType SNP

	java -Xmx96G -jar $GATK \
		-T SelectVariants \
		-R $refGenome \
		-V chr${d}.vcf \
		-o chr${d}_indels.vcf \
		-selectType INDEL
}

hardfiltering_snvs () {
	d=$1
	java -Xmx96G -jar $GATK \
		-T VariantFiltration \
		-R $refGenome \
		-V chr${d}_snvs.vcf \
		--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
		--filterName "PAD_snv_filter" \
		-o chr${d}_snvs_filtered.vcf
}

hardfiltering_indels () {
	d=$1
	java -Xmx96G -jar $GATK \
		-T VariantFiltration \
		-R $refGenome \
		-V chr${d}_indels.vcf \
		--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
		--filterName "PAD_indel_filter" \
		-o chr${d}_indels_filtered.vcf
}

concat_vcf () {

	vcflist=`find ./chr* -name "chr*_snvs_filtered.vcf" -exec echo "-V " {} \;`
	java -jar $GATK -T CombineVariants \
		-R $refGenome \
		$vcflist \
		-o $OUT_DIR/PAD_latest.SNVs.vcf \
		-genotypeMergeOptions UNIQUIFY

	vcflist2=`find ./chr* -name "chr*_indels_filtered.vcf" -exec echo "-V " {} \;`
	java -jar $GATK -T CombineVariants \
		-R $refGenome \
		$vcflist2 \
		-o $OUT_DIR/PAD_latest.indels.vcf \
		-genotypeMergeOptions UNIQUIFY

}

concat_keyue () {

	vcflist=`find ./chr* -name "chr*_snvs_filtered.vcf" -exec echo "-V " {} \;`
	java -cp $GATK org.broadinstitute.gatk.tools.CatVariants \
		-R $refGenome \
		$vcflist \
		-out $OUT_DIR/PAD_latest.SNVs.vcf \
		-assumeSorted

	vcflist2=`find ./chr* -name "chr*_indels_filtered.vcf" -exec echo "-V " {} \;`
	java -cp $GATK org.broadinstitute.gatk.tools.CatVariants \
		-R $refGenome \
		$vcflist2 \
		-out $OUT_DIR/PAD_latest.indels.vcf \
		-assumeSorted


}

# # # main # # # 
#

source $HOME/.bash_profile

cd /data5/dom/rs-cardio/kullo/emerge2/users/m003388/PAD_seq/combinedvcf 

#chr=$1

#cd chr${chr}

#combineGVCFs $chr

#genotypeGVCFs $chr

#selectVariants $chr 

#hardfiltering_snvs $chr 

#hardfiltering_indels $chr

concat_keyue
