#! /bin/sh

# Variant Calling
# Institute: Mayo clinic 

# reference genome (b37), dbsnp
gatk_db=/data5/dom/rs-cardio/kullo/emerge2/dataBase/public/resource_bundle
refGenome=$gatk_db/human_g1k_v37_decoy.fasta
dbsnp=$gatk_db/dbsnp_138.b37.vcf
known_indelsvcf=$gatk_db/1000G_phase1.indels.b37.vcf
known_millsvcf=$gatk_db/Mills_and_1000G_gold_standard.indels.b37.vcf

# programs: bwa (ver: 0.7.13), GATK (ver: 3.5), and picard (2.8.0)
software_dir=/data5/dom/rs-cardio/kullo/emerge2/softWare
bwa=$software_dir/bwa-0.7.12/bwa
GATK=$software_dir/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
picard=$software_dir/picard/build/libs/picard.jar

# original fastq files
originaldir=/data5/dom/rs-cardio/kullo/emerge2/dataBase/sequencing/PADseq/data/bamFiles
keyuedir=/data5/dom/rs-cardio/kullo/emerge2/users/m003388/PAD_seq/haplotypeCaller/$1

# temporary folder
TMP_DIR=/local2/tmp/m164988


# # # # # # # # # # # # #
# # # # functions # # # #
# # # # # # # # # # # # #

markIlluminaAdapters () {
	sample_id=$1
	# mark adapter sequences using MarkIlluminaAdapters
	java -Xmx16G -jar $picard MarkIlluminaAdapters \
		I=${sample_id}.sorted.bam \
		O=${sample_id}_markilluminaadapters.bam \
		M=${sample_id}_markilluminaadapters_metrics.txt \
		TMP_DIR=$TMP_DIR
}

Samtofastq () {
	sample_id=$1
	# Convert BAM to FASTQ and discount adapter sequences using SamToFastq
	java -Xmx16G -jar $picard SamToFastq \
		I=${sample_id}_markilluminaadapters.bam \
		FASTQ=${sample_id}_samtofastq_interleaved.fq \
		CLIPPING_ATTRIBUTE=XT \
		CLIPPING_ACTION=2 \
		INTERLEAVE=true \
		NON_PF=true \
		TMP_DIR=$TMP_DIR
}

bwa_mem () {
	sample_id=$1
	# alignment using bwa-mem
	$bwa mem -M -t 4 -p $refGenome ${sample_id}_samtofastq_interleaved.fq > ${sample_id}_bwa_mem.sam

	# -M to flag shorter split hits as secondary. This is optional for Picard compatibility as MarkDuplicates can directly process BWA's alignment, whether or not the alignment marks secondary hits. However, if we want MergeBamAlignment to reassign proper pair alignments, to generate data comparable to that produced by the Broad Genomics Platform, then we must mark secondary alignments.
	# -p to indicate the given file contains interleaved paired reads.
	# -t followed by a number for the number of processor threads to use concurrently. Here we use seven threads which is one less than the total threads available on my Mac laptap. Check your server or system's total number of threads with the following command provided by KateN.
}

mergeBamAlignment () {
	sample_id=$1
	# merge with uBAM
	java -Xmx16G -jar $picard MergeBamAlignment \
		R=$refGenome \
		UNMAPPED_BAM=${sample_id}.sorted.bam \
		ALIGNED_BAM=${sample_id}_bwa_mem.sam \
		O=${sample_id}_mergebamalignment.bam \
		CREATE_INDEX=true \
		ADD_MATE_CIGAR=true \
		CLIP_ADAPTERS=false \
		CLIP_OVERLAPPING_READS=true \
		INCLUDE_SECONDARY_ALIGNMENTS=true \
		MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
		ATTRIBUTES_TO_RETAIN=XS \
		TMP_DIR=$TMP_DIR
}

markDuplicates () {
	sample_id=$1
	# mark duplicates and remove 
	java -Xmx16G -jar $picard MarkDuplicates \
		INPUT=${sample_id}_mergebamalignment.bam \
		OUTPUT=${sample_id}_markduplicates.bam \
		METRICS_FILE=${sample_id}_markduplicates_metrics.txt \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
		REMOVE_DUPLICATES=true \
		CREATE_INDEX=true \
		TMP_DIR=$TMP_DIR
}

realignTargetCreater () {
	# create indels target region
	sample_id=$1
	java -Xmx16G -jar $GATK \
		-T RealignerTargetCreator \
		-R $refGenome \
		-known $known_indelsvcf \
		-known $known_millsvcf \
		-I ${sample_id}_markduplicates.bam \
		-o ${sample_id}_realignertargetcreator.intervals
}


indelRealigner () {
	# indels realignment
	sample_id=$1
	java -Xmx16G -Djava.io.tmpdir=$TMP_DIR -jar $GATK -T IndelRealigner \
		-R $refGenome \
		-targetIntervals ${sample_id}_realignertargetcreator.intervals \
		-known $known_indelsvcf \
		-known $known_millsvcf \
		-I ${sample_id}_markduplicates.bam \
		-o ${sample_id}_indelrealigner.bam
}

run_BQSR () {
	# Recalibrate base quality scores
	sample_id=$1

	# Analyze patterns of covariation in the sequence dataset
	java -Xmx16G -jar $GATK -T BaseRecalibrator \
		-R $refGenome \
		-I ${sample_id}_indelrealigner.bam \
		-knownSites $dbsnp \
		-knownSites $known_indelsvcf \
		-knownSites $known_millsvcf \
		-o recal_data.table
 
	# Do a second pass to analyze covariation remaining after recalibration
	java -Xmx16G -jar $GATK -T BaseRecalibrator \
		-R $refGenome \
		-I ${sample_id}_indelrealigner.bam \
		-knownSites $dbsnp \
		-knownSites $known_indelsvcf \
		-knownSites $known_millsvcf \
		-BQSR recal_data.table \
		-o post_recal_data.table 
#	
#	# Generate before/after plot
#	java -Xmx64G -jar $GATK -T AnalyzeCovariates \
#		-R $refGenome \
#		-before recal_data.table \
#		-after post_recal_data.table \
#		-plots recalibration_plots.pdf

	# Apply the recalibration to your sequence data 
	java -Xmx16G -jar $GATK -T PrintReads \
		-R $refGenome \
		-I ${sample_id}_indelrealigner.bam \
		-BQSR recal_data.table \
		-o ${sample_id}_recal_reads.bam 
}

haplotyperCaller () {
	sample_id=$1
	java -Xmx64G -jar $GATK -T HaplotypeCaller \
		-R $refGenome \
		-ERC GVCF \
		-I $keyuedir/${sample_id}_recal_reads.bam \
		--genotyping_mode DISCOVERY \
		--variant_index_type LINEAR \
		--variant_index_parameter 128000 \
		-L ../kullo_rsng_r216_nimblegen_469kb.list \
		-stand_call_conf 10 \
		-o ${sample_id}_raw_targets_variants.g.vcf
}

unzipBam () {
	sample_id=$1
	gunzip -c $originaldir/${sample_id}.bam.gz > ./${sample_id}_recal_reads.bam
	gunzip -c $originaldir/${sample_id}.bai.gz > ./${sample_id}_recal_reads.bai
}


select_variants () {
	sample_id=$1
	java -Xmx64G -jar $GATK -T SelectVariants \
		-R $refGenome \
		-L ../kullo_rsng_r216_nimblegen_469kb.list \
		-V ${sample_id}_raw_variants.g.vcf \
		-o ${sample_id}_raw_variants.targeted.g.vcf
}

# # # # # # # # #  # # #
# # # main program # # #
# # # # # # # # #  # # #

sample_id=$1

if ! [ -d $sample_id ]; then mkdir $sample_id; fi; cd $sample_id

# # 1. Convert Fastq to uBAM and add read group information using FastqToSam 
#fastq2uBAM $sample_id
#unzip_sort_Bam $sample_id

## # 2. Mark adapter sequences using MarkIlluminaAdapters
#markIlluminaAdapters $sample_id

## # 3. Convert BAM to FASTQ and discount adapter sequences using SamToFastq
#Samtofastq $sample_id

## # 4. align with BWA_MEM and merge with uBAM
#bwa_mem $sample_id
#mergeBamAlignment $sample_id

## # 5. mark duplicates
#unzipBam $sample_id
#markDuplicates $sample_id

## # 6. realignment
#realignTargetCreater $sample_id
#indelRealigner $sample_id

## # 7. Recalibrate base quality scores in order to correct sequencing errors and other experimental artifacts
#run_BQSR $sample_id

## # 8. variant callng
haplotyperCaller $sample_id

#select_variants $sample_id
