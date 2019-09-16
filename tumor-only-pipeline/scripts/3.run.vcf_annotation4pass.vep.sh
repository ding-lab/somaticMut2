
# Hua Sun
# 9/12/2019

# VEP annotation
# the parameter followed somaticwrapper pipeline


# sh run.sh -i <vcf>
# sh run.sh -i <vcf> -a <GRCh38> -g <GRCh38.genome.fa> -d <vep_cache_dir>
# output same as input file & the out name as *.vep_anno



# Set in MGI-Server
#export PERL_PATH=/gsc
#export PATH=$PERL_PATH/bin:$PATH
#export PERL_BIN=$PERL_PATH/bin/perl
#export PERL5LIB=$PERL_PATH/lib/perl5/5.8.7/:$PERL5LIB


while getopts "c:i:a:g:d:" opt; do
  case $opt in
  	c)
  		CONFIG=$OPTARG
      ;;
    i)
    	VCF=$OPTARG
      ;;
    a)
    	ASSEMBLY=$OPTARG
      ;;    	
    g)
    	GENOME=$OPTARG
    	;;
    d)
    	VEP_CACHE_DIR=$OPTARG
    	;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


source $CONFIG


# format VCF for VEP
grep ^# ${VCF} > ${VCF}.header

# extract 'PASS' variant sites
# For Mutect2 DP >= 5
grep -v ^# ${VCF} | awk -F['\t'] '$7=="PASS"' | perl -ne '$DP=$1 if /DP=(\d+)/; print if $DP>=5' > ${VCF}.fmt4vep
# with 'chr' or without 'chr' both okay


# if use conda version, no need to perl 
perl $VEP --buffer_size 300 --offline --cache \
		--dir ${VEP_CACHE_DIR} \
		--assembly ${ASSEMBLY} \
		--fasta ${GENOME} \
		--fork 8 --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol \
		--numbers --domains --canonical --protein --biotype --uniprot --tsl --pubmed \
		--variant_class --shift_hgvs 1 --check_existing --total_length --allele_number \
		--no_escape --xref_refseq --failed 1 --minimal --flag_pick_allele --force_overwrite \
		--pick_order canonical,tsl,biotype,rank,ccds,length \
		--format vcf --vcf -i ${VCF}.fmt4vep -o ${VCF}.fmt4vep.vep_anno
		
fileName=`echo ${VCF} | sed 's/\.vcf$//'`
mv ${VCF}.fmt4vep.vep_anno ${fileName}.vep.vcf

# remove
rm -f ${VCF}.header ${VCF}.fmt4vep


