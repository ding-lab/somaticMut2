=head1

  
  Hua Sun
  v0.1 9/5/2019
  
  
  ToDo
  	Extract filtered vcf from Mutect2 tumor-only VCF file
  	Format VCF to TSV (Form2 format)
  
  
  Usage
    perl vcf_filter_and_makeTSV.m2snpeff.pl mutect2.vcf > filtered.pass.tsv
    perl vcf_filter_and_makeTSV.m2snpeff.pl -proteinCoding mutect2.vcf > protein_coding.filtered.pass_anno.tsv
    
    -minVaf   [float]    Minimum VAF (default 0.05)
    -minDepth [int]      Minimum Depth (default 8)  
    -minMut   [int]      Minimum Mutation (default 3)
    
    -proteinCoding       only call protein coding region
                         indel, snv, stop_gain/lost, splice
    


  Input
  	mutect2.snpEff.tsv
  
  Output
    Chr	Start_Pos	Ref	Alt	Gene	Variant	VariantType	Sample	T_Depth	T_RefCount	T_AltCount	T_VAF	Genotype	Info	HGVSp	Classification
    chr14	24131583	T	G	FITM1	p.V7G	Missense_Mutation	H_MT-389-D34A_1674	69	62	7	0.101449275	0/1	SNP;HET	p.Val7Gly	missense_variant    

 
=cut


use strict;
use Getopt::Long;

my $min_vaf = 0.05;
my $min_depth = 8;
my $min_mut = 3;
my $protein_coding;
my $Help;

GetOptions(
  "minVaf:f" => \$min_vaf,
  "minDepth:i" => \$min_depth,
  "minMut:i" => \$min_mut,
  "proteinCoding" => \$protein_coding,
  "help" => \$Help
);


die `pod2text $0` if (@ARGV==0 || $Help);

# chr1  13390963  .   C   G   .   PASS  CONTQ=93;DP=75;ECNT=1;GERMQ=14;MBQ=20,20;MFRL=162,136;MMQ=60,48;MPOS=33;POPAF=7.30;ROQ=93;SEQQ=31;STRANDQ=64;TLOD=7.99  GT:AD:AF:DP:F1R2:F2R1:SB        0/1:63,8:0.121:71:35,4:26,4:25,38,4,4

my $file = shift;
my $header = `grep '#CHROM' $file`;

my $sampleName = `grep '#CHROM' $file | cut -f 10`;
chomp $sampleName;

my @data = `grep -v '^#' $file | awk -F['\t'] '\$7=="PASS"'`;

# header
print "Chr	Start_Pos	End_Pos	Ref	Alt	Gene	Variant	VariantType	Sample	T_Depth	T_RefCount	T_AltCount	T_VAF	N_Depth	N_RefCount	N_AltCount	N_VAF	Genotype	Info	HGVSp	Classification\n";



# main
my ($chr, $start_pos, $end_pos);
my ($ref, $alt);
my ($GT, $AD, $DP, $VAF);
my ($ref_count, $alt_count);
my ($gene, $mut_aa3, $mut_aa1, $mut_type, $mut_type_new, $info);


my $anno_info;
my $pos;

foreach (@data) {
	
	my @arr = split("\t");
	
	$chr = $arr[0];
	$start_pos = $arr[1];
	$ref = $arr[3];
	$alt = $arr[4];
	
	# cal. end_pos
	if (length($ref) == length($alt)){
	    # SNP
	    if (length($ref) == 1){
	        $end_pos = $start_pos;
	    } else {
	    # MNP
	        $end_pos = $start_pos + length($ref) - 1;
	    }
		
	} elsif (length($ref) > length($alt)) {
	    $end_pos = $start_pos + length($ref) - 1;
	    
	} else {
			$end_pos = $start_pos + 1;
			
	}
	
	
	# sample INFO
	# GT:AD:AF:DP:F1R2:F2R1:SB
	# 0/1:63,8:0.121:71:35,4:26,4:25,38,4,4
	my @arr_sampleInfo = split(":", $arr[9]);
	
	$GT = $arr_sampleInfo[0];
	$AD = $arr_sampleInfo[1];
	($ref_count, $alt_count) = split(",", $AD);
	$DP = $arr_sampleInfo[3];
	$VAF = $alt_count/$DP;
	
	$anno_info = $arr[7];
	
  ## Filter by DP, ALT, VAF
  if ($DP >= $min_depth && $alt_count >= $min_mut && $VAF >= $min_vaf){
  
       # extract anno info from snpEff
      ($gene, $mut_aa3, $mut_aa1, $mut_type, $mut_type_new, $info) = ExtractInfo_from_snpEffAnno($anno_info);

  	  ## Coding mutation
      if ( $protein_coding ){
      
          next unless ($mut_type_new =~ /Frame/ || $mut_type_new eq 'Splice_Site' || $mut_type_new =~ /_Mutation/ || $mut_type_new eq 'Silent')
		
          # http://snpeff.sourceforge.net/SnpEff_manual.html
          #my $var_type = (split(/\|/, $anno_info))[1];
          #next unless ($var_type =~ /inframe_/ || $var_type =~ /frameshift_/ || $var_type =~ /splice_acceptor/ || $var_type =~ /splice_donor/ || $var_type =~ /exon_loss/ || $var_type =~ /stop_/ || $var_type =~ /missense_/ || $var_type =~ /coding_sequence/ || $var_type =~ /rare_amino_acid/ || $var_type =~ /synonymous/);
      }
      
      print "$chr\t$start_pos\t$end_pos\t$ref\t$alt\t$gene\t$mut_aa1\t$mut_type_new\t$sampleName\t$DP\t$ref_count\t$alt_count\t$VAF\t\.\t\.\t\.\t\.\t$GT\t$info\t$mut_aa3\t$mut_type\n";
  }
  
}


exit;




## ExtractInfo_from_snpEffAnno
sub ExtractInfo_from_snpEffAnno
{
	my $str_anno = shift;
	
	my @arr = split(/\|/, $str_anno);
	
	my ($Hugo_Symbol, $HGVSp, $Variant_Classification);

	# INFO
	my $Info = $arr[0];
	my $hom_het = '';
	my $var_type = '';
	if ($Info=~/\;(\w+)\;VARTYPE=(\w+)\;/){
		$hom_het = $1;
		$var_type = $2;
	}

	
	$Variant_Classification = $arr[1];
	$Hugo_Symbol = $arr[3];
	$HGVSp = $arr[10];
	$HGVSp = '.' if ($HGVSp eq '');
	
	my $HGVSp_Short = AminoAcid_3to1( $HGVSp );
	my $Variant_Classification_new = GetVariantClassification( $Variant_Classification, $var_type, $HGVSp );

	
	my $func = '';
	$func = 'LOF' if ($str_anno =~ /\;LOF\=/);
	if ($str_anno =~ /\:Phosphotyrosine\|/){
		if ($func ne ''){
			$func .= ';Phosphotyrosine';
		} else {
			$func = 'Phosphotyrosine';
		}
	}
	
	if ($func eq ''){
		$Info = "$var_type\;$hom_het";
	} else {
		$Info = "$var_type\;$hom_het\;$func";
	}
	
	
	return ($Hugo_Symbol, $HGVSp, $HGVSp_Short, $Variant_Classification, $Variant_Classification_new, $Info);
}



## Converts Sequence Ontology variant types to MAF variant classifications
## The original code from https://github.com/mskcc/vcf2maf/blob/master/vcf2maf.pl
## Re-writed code to fit the SnpEff Annotation
sub GetVariantClassification
{
		#my ( $effect, $var_type, $inframe) = @_;
		my ($effect, $var_type, $HGVSp) = @_;
		
		# use to SnpEff coding region
    return "In_Frame_Ins" if( $effect =~ /inframe_insertion/ and $HGVSp ne '.' );
    return "In_Frame_Del" if( $effect =~ /inframe_deletion/ and $HGVSp ne '.' );
    return "Frame_Shift_Del" if( $effect =~ /frameshift_variant/ and $var_type eq 'DEL' and $HGVSp ne '.' );
    return "Frame_Shift_Ins" if( $effect =~ /frameshift_variant/ and $var_type eq 'INS' and $HGVSp ne '.' );
    return "Splice_Site" if( ($effect eq 'splice_acceptor_variant' || $effect eq 'splice_donor_variant' || $effect eq 'exon_loss_variant') and $HGVSp eq '.');
    return "Nonsense_Mutation" if( $effect =~ /stop_gained/ and $HGVSp ne '.' );   
    return "Nonstop_Mutation" if( $effect =~ /stop_lost/ and $HGVSp ne '.' );
    return "Missense_Mutation" if( $effect =~ /missense_variant/ || $effect eq 'coding_sequence_variant' || $effect eq 'rare_amino_acid_variant' );
    return "Silent" if( $effect eq 'synonymous_variant' || $effect eq 'stop_retained_variant' );
    return "Splice_Region" if( $effect =~ /splice_region_variant/ and $HGVSp eq '.');
    return "Translation_Start_Site" if( $effect eq 'initiator_codon_variant' || $effect eq 'start_lost' );
    return "Intron" if( $effect =~ 'intron_variant' || $effect eq 'initiator_codon_variant' || $effect eq 'intragenic_variant');
    return "IGR" if( $effect =~ /intergenic_/ || $effect eq 'TF_binding_site_variant' || $effect eq 'structural_interaction_variant' || $effect eq 'protein_protein_contact' );
    return "RNA" if( $effect =~ /non_coding/ );
    return "5'UTR" if( $effect =~ /5_prime_UTR/ );
    return "3'UTR" if( $effect =~ /3_prime_UTR/ );
    return "5'Flank" if( $effect eq 'upstream_gene_variant' );
    return "3'Flank" if( $effect eq 'downstream_gene_variant' );
    return "Sequence_feature" if( $effect eq 'sequence_feature' );

    
    # original for VEP
    #return "Splice_Site" if( $effect =~ /^(splice_acceptor_variant|splice_donor_variant|transcript_ablation|exon_loss_variant)$/ );
    #return "Nonsense_Mutation" if( $effect eq 'stop_gained' );
    #return "Frame_Shift_Del" if(( $effect eq 'frameshift_variant' or ( $effect eq 'protein_altering_variant' and !$inframe )) and $var_type eq 'DEL' );
    #return "Frame_Shift_Ins" if(( $effect eq 'frameshift_variant' or ( $effect eq 'protein_altering_variant' and !$inframe )) and $var_type eq 'INS' );
    #return "Nonstop_Mutation" if( $effect eq 'stop_lost' );
    #return "Translation_Start_Site" if( $effect =~ /^(initiator_codon_variant|start_lost)$/ );
    #return "In_Frame_Ins" if( $effect =~ /^(inframe_insertion|disruptive_inframe_insertion)$/ or ( $effect eq 'protein_altering_variant' and $inframe and $var_type eq 'INS' ));
    #return "In_Frame_Del" if( $effect =~ /^(inframe_deletion|disruptive_inframe_deletion)$/ or ( $effect eq 'protein_altering_variant' and $inframe and $var_type eq 'DEL' ));
    #return "Missense_Mutation" if( $effect =~ /^(missense_variant|coding_sequence_variant|conservative_missense_variant|rare_amino_acid_variant)$/ );
    #return "Intron" if ( $effect =~ /^(transcript_amplification|intron_variant|INTRAGENIC|intragenic_variant)$/ );
    #return "Splice_Region" if( $effect eq 'splice_region_variant' );
    #return "Silent" if( $effect =~ /^(incomplete_terminal_codon_variant|synonymous_variant|stop_retained_variant|NMD_transcript_variant)$/ );
    #return "RNA" if( $effect =~ /^(mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|nc_transcript_variant)$/ );
    #return "5'UTR" if( $effect =~ /^(5_prime_UTR_variant|5_prime_UTR_premature_start_codon_gain_variant)$/ );
    #return "3'UTR" if( $effect eq '3_prime_UTR_variant' );
    #return "IGR" if( $effect =~ /^(TF_binding_site_variant|regulatory_region_variant|regulatory_region|intergenic_variant|intergenic_region)$/ );
    #return "5'Flank" if( $effect eq 'upstream_gene_variant' );
    #return "3'Flank" if ( $effect eq 'downstream_gene_variant' );


    # Annotate everything else simply as a targeted region
    # TFBS_ablation, TFBS_amplification,regulatory_region_ablation, regulatory_region_amplification,
    # feature_elongation, feature_truncation
    return "Targeted_Region";
}



## Convert 3-letter amino-acid to 1-letter
## ref: https://github.com/mskcc/vcf2maf/blob/master/vcf2maf.pl
sub AminoAcid_3to1
{
    my $aa3 = shift;
    
    # Hash to convert 3-letter amino-acid codes to their 1-letter codes
    my %aa3to1 = qw( Ala A Arg R Asn N Asp D Asx B Cys C Glu E Gln Q Glx Z Gly G His H Ile I Leu L
                     Lys K Met M Phe F Pro P Ser S Thr T Trp W Tyr Y Val V Xxx X Ter * );
    
   my $aa1 = $aa3;
     
	 # Create a shorter HGVS protein format using 1-letter codes
	  while( $aa1 and my ( $find, $replace ) = each %aa3to1 ) {
        eval "\$aa1 =~ s{$find}{$replace}g";
    }
    
    return $aa1;
}

