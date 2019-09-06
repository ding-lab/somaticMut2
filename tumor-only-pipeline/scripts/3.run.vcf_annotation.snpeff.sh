# 8/18/2019
# annotation for vcf
# min memory 16GB - 5 min

VCF=$1


JAVA=/gscmnt/gc2737/ding/hsun/software/jre1.8.0_152/bin/java
SNPEFF=/gscmnt/gc3021/dinglab/hsun/software/snpEff/snpEff.jar
SNPSIFT=/gscmnt/gc3021/dinglab/hsun/software/snpEff/SnpSift.jar


# mutation annotation
# PASS + HOM/HET + LoF
cat $VCF | $JAVA -Xmx16G -jar ${SNPSIFT} filter "(FILTER = 'PASS')" | $JAVA -Xmx16G -jar ${SNPSIFT} varType - | $JAVA -Xmx16G -jar $SNPEFF GRCh38.86 - > ${VCF}.pass.anno


# anno all mutation sites
#$JAVA -Xmx16G -jar $SNPEFF -v GRCh38.86 $VCF > ${VCF}.anno

# \;ANN=C|<annotation>|

