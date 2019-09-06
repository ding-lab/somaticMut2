
# Hua Sun
# 8/13/2019

# USAGE
# sh run.sh -i <vcf> -d <dbSNP_noCOSMIC>
# *.vcf (not *.vcf.gz)

# http://snpeff.sourceforge.net/SnpSift.html#filter

JAVA=/gscmnt/gc2737/ding/hsun/software/jre1.8.0_152/bin/java
SNPSIFT=/gscmnt/gc3021/dinglab/hsun/software/snpEff/SnpSift.jar
dbSNP_noCOSMIC=/gscmnt/gc3021/dinglab/hsun/Database/dbSNP/All_20180418.b151.filtered.CosmicCoding_NonCoding_v89.GRCh38.vcf.gz

# set memory 8GB

while getopts "i:d:" opt; do
  case $opt in
    i)
    	vcf=$OPTARG
      ;;
    d)
    	dbSNP_noCOSMIC=$OPTARG
    	;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


if [ ! -f $dbSNP_noCOSMIC ]; then
	echo "[ERROR] No dbSNP file or direction ..." >&2
	exit 1
fi

echo "[INFO] Filter dbSNP location from vcf ..."  >&2

file=${vcf%.vcf}
OUT=$file.rem_dbSNP_noCOSMIC.vcf


# pre-filter for saving running time of SnpSift
grep -v '^#' $vcf | awk -F['\t'] '$7=="PASS"' > $vcf.pass.tmp
grep '^#' $vcf > $vcf.header.tmp
cat $vcf.header.tmp $vcf.pass.tmp > $vcf.pass.h.tmp

vcf_input=$vcf.pass.h.tmp


# filter dbSNP
$JAVA -Xmx8G -jar $SNPSIFT annotate -id $dbSNP_noCOSMIC $vcf_input | $JAVA -Xmx8G -jar $SNPSIFT filter -n "(exists ID) & (ID =~ 'rs' )" > $OUT


# remove tempfile
rm -f $vcf.pass.tmp $vcf.header.tmp $vcf.pass.h.tmp

