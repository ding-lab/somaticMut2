
# Hua Sun
# Version 1.1
# 8/11/2019


# Mutect 2 - pipeline
# https://software.broadinstitute.org/gatk/documentation/article?id=24057


# USAGE
# sh run.sh -n <sampleName> -b <bam> -o <outdir>

# Memory 4Gb, running 9h in MGI per sample
# set 8Gb for final step


# Setting

GATK=/gscmnt/gc3021/dinglab/hsun/software/gatk-4.1.2.0/gatk

SAMTOOLS=/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/samtools

JAVA=/gscmnt/gc2737/ding/hsun/software/jre1.8.0_152/bin/java
PICARD=/gscmnt/gc2737/ding/hsun/software/picard.jar


# GECODE reference
ref_fa=/gscmnt/gc2737/ding/hsun/data/human_genome/gencode_GRCh38_v29/genome/GRCh38.primary_assembly.genome.fa
# GDC reference
#ref_fa=/gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data/image.data/A_Reference/GRCh38.d1.vd1.fa

# gnomad
gnomad_vcf=/gscmnt/gc3021/dinglab/hsun/Database/gatk/mutect2_gatk-best-practices.broadIns/af-only-gnomad.hg38.vcf.gz

# pon normal
panel_of_normals_vcf=/gscmnt/gc3021/dinglab/hsun/Database/gatk/GDC-gatk4_panel_of_normals/6c4c4a48-3589-4fc0-b1fd-ce56e88c06e4/gatk4_mutect2_4136_pon.hg38.vcf.gz

# common_biallelic
common_biallelic=/gscmnt/gc3021/dinglab/hsun/Database/gatk/mutect2_gatk-best-practices.broadIns/af-only-gnomad.hg38.common_biallelic.chr1-22XY.vcf


name=''
outdir=`pwd`/mutect2_result


while getopts "n:b:r:o:" opt; do
  case $opt in
    n)
    	name=$OPTARG
      ;;  
    b)
    	tumor_bam=$OPTARG
      ;;
    r)
    	ref_fa=$OPTARG
    	;;
    o)
      outdir=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done



# check input
if [[ $name == '' ]];then
  echo "[ERROR] Please input -n sampleName ..." >&2
  exit
fi

if [ ! -f "$tumor_bam" ];then
  echo "[ERROR] The $tumor_bam is not existing ..." >&2
  exit
fi



# make outdir folder
mkdir -p $outdir

OUT=$outdir/$name
mkdir -p $OUT


##==================== Make softlink for bam ====================##
# Create softlink bam
if [ ! -f "$OUT/$name.bam" ];then
  ln -s $tumor_bam $OUT/$name.bam
fi

if [ -f "$tumor_bam.bai" ];then
  ln -s $tumor_bam.bai $OUT/$name.bam.bai
else
  $SAMTOOLS index $OUT/$name.bam
fi





##==================== 1. Make a unfilter VCF ====================##
## Make unfilter file
# *f1r2.tar.gz  *unfiltered.vcf/idx/status
for chromosome in {'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'}; do
  
  $GATK --java-options "-Xmx4g" Mutect2 \
     -I $OUT/$name.bam \
     -R ${ref_fa} \
     -L ${chromosome} \
     --germline-resource ${gnomad_vcf} \
     -pon ${panel_of_normals_vcf} \
     --f1r2-tar-gz ${OUT}/${chromosome}-f1r2.tar.gz \
     -O ${OUT}/${chromosome}-unfiltered.vcf
 
done



## Merge unfiltered-vcf
all_unfiltered_input=`for chrom in {1..22}; do printf -- "I=${OUT}/chr${chrom}-unfiltered.vcf "; done; for chrom in {'X','Y'}; do printf -- "I=${OUT}/chr${chrom}-unfiltered.vcf "; done`

$JAVA -Xmx4G -jar $PICARD GatherVcfs \
    $all_unfiltered_input \
    O=${OUT}/merged-unfiltered.vcf


# Merged vcf.stats
all_unfiltered_stats_input=`for chrom in {1..22}; do printf -- "-stats ${OUT}/chr${chrom}-unfiltered.vcf.stats "; done; for chrom in {'X','Y'}; do printf -- "-stats ${OUT}/chr${chrom}-unfiltered.vcf.stats "; done`

$GATK MergeMutectStats \
    $all_unfiltered_stats_input \
    -O ${OUT}/merged-unfiltered.vcf.stats





##==================== 2. Make contamination table ====================##
## Make read-orientation-model
# make 'read-orientation-model.tar.gz'
all_f1r2_input=`for chrom in {1..22}; do printf -- "-I ${OUT}/chr${chrom}-f1r2.tar.gz "; done; for chrom in {'X','Y'}; do printf -- "-I ${OUT}/chr${chrom}-f1r2.tar.gz "; done`
	
# it must write to like this
$GATK --java-options "-Xmx4g" LearnReadOrientationModel ${all_f1r2_input} -O ${OUT}/read-orientation-model.tar.gz


## Make 'getpileupsummaries.table'
$GATK --java-options "-Xmx4g" GetPileupSummaries \
    -I $OUT/$name.bam \
    -V ${common_biallelic} \
    -L ${common_biallelic} \
    -O ${OUT}/getpileupsummaries.table


## Make 'contamination.table' & 'segments.table'
$GATK --java-options "-Xmx4g" CalculateContamination \
    -I ${OUT}/getpileupsummaries.table \
    --tumor-segmentation ${OUT}/segments.table \
    -O ${OUT}/contamination.table



# 8/25/2019
# NOTE: the FilterMutectCalls will not run when the contamination.table vale is 'NaN'
# Solution is change the value to 0 0
perl -i -pe 's/\tNaN\t1\.0/\t0\t0/ if /\tNaN\t1\.0/' ${OUT}/contamination.table


##==================== 3. Make a final filter VCF ====================##
# Make a finial filtered vcf
$GATK --java-options "-Xmx8g" FilterMutectCalls \
    -V ${OUT}/merged-unfiltered.vcf \
    -R ${ref_fa} \
    --tumor-segmentation ${OUT}/segments.table \
    --contamination-table ${OUT}/contamination.table \
    --ob-priors ${OUT}/read-orientation-model.tar.gz \
     -O ${OUT}/filtered.vcf




##==================== remove temporary file ====================##
rm -f ${OUT}/chr*




