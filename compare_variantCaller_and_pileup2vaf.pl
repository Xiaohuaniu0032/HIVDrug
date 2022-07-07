use strict;
use warnings;
use Data::Dumper;

my ($allele_file,$pileup_file) = @ARGV;

# my $freq_cutoff = 0.04;
# /data/fulongfei/analysis/hiv/re_analysis_new_ref_with_hs/IonXpress_002_rawlib/variantCaller

my %allele_var;
my %hs_snv_by_order;

my %hs_pos;
open IN, "$allele_file" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $pos = $arr[1];
	my $ref = $arr[2];
	my $alt = $arr[3];
	my $freq = $arr[6];

	if (/Hotspot/){
		$hs_pos{$arr[1]}{$arr[2]}{$arr[3]} = 1; # pos/ref/alt
	}
}
close IN;

open IN, "$allele_file" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $pos = $arr[1];
	my $ref = $arr[2];
	my $alt = $arr[3];
	my $freq = $arr[6];

	#next if ($arr[4] eq "Absent" || $arr[4] eq "No Call");

	# HS SNP
	#if ($arr[10] eq "Hotspot"){
	#	my $var = "$arr[1]\t$arr[2]\t$arr[3]"; # pos/ref/alt
	#	push @{$hs_snv_by_order{$arr[1]}}, "$arr[1]\t$arr[2]\t$arr[3]"; # {pos=>[pos/ref/alt,...,]}
	#}

	# alleles.xls文件中，有些位点会有重复
	# HIV-1   2806    T       C       Heterozygous    -       38.5    5864.9  -       SNP     Novel
	# HIV-1   2806    TAA     CAG     Heterozygous    -       60.9    5864.9  -       MNP     Novel

	my $val = "$arr[4]\t$freq\t$arr[10]"; # Heterozygous/6.7/Hotspot
	
	my $type = $arr[9]; # MNP/SNP
	if ($type eq "SNP"){
		my $var = "$arr[1]\t$arr[2]\t$arr[3]"; # pos/ref/alt
		push @{$allele_var{$var}}, $val;
	}elsif ($type eq "MNP"){
		my $len = length($ref);
		my @ref_arr = split //, $ref;
		my @alt_arr = split //, $alt;
		for my $i (0..$len-1){
			my $ref_base = $ref_arr[$i];
			my $alt_base = $alt_arr[$i];
			if ($ref_base ne $alt_base){
				my $rel_pos = $pos + $i;
				my $var = "$rel_pos\t$ref_base\t$alt_base"; # 如果热点SNP未检出但Novel检出则会显示两行
				push @{$allele_var{$var}}, $val;
			}
		}
	}else{
		next; # do not consider INDEL/COMPLEX
	}
}
close IN;
#print(Dumper(\%allele_var));


my %pileup_var;
# Sample  Chr     Pos     Ref     Alt     AltNum  Depth   AltAlleleFrequency
open IN, "$pileup_file" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $ref = $arr[3];
	my $alt = $arr[4];
	
	my $ref_len = length($ref);
	my $alt_len = length($alt);
	next if ($ref_len != $alt_len); # skip indel
	
	my $var = "$arr[2]\t$ref\t$alt"; # pos/ref/alt
	$pileup_var{$var} = $arr[-1] * 100;
}
close IN;


print "Chrom\tPos\tRef\tAlt\tAllele_Call\tTVC_Freq\tTVC_Note\tPileup_Freq\tif_PASS\tif_Report\n";

foreach my $pos (sort {$a <=> $b} keys %hs_pos){
	# $hs_pos{$arr[1]}{$arr[2]}{$arr[3]} = 1; # pos/ref/alt
	my @ref = keys %{$hs_pos{$pos}};
	for my $ref (@ref){
		my @alt = keys %{$hs_pos{$pos}{$ref}};
		for my $alt (@alt){
			my $var = "$pos\t$ref\t$alt";

			# pileup info
			my $pileup_freq;
			if (exists $pileup_var{$var}){
				$pileup_freq = $pileup_var{$var};
			}else{
				$pileup_freq = 0;
			}

			# tvc info
			my @tvc_info = @{$allele_var{$var}};
			for my $info (@tvc_info){
				# Heterozygous/6.7/Hotspot
				my $tvc_freq = (split /\t/, $info)[1];

			
				my $report_cutoff = 3; # 3%
				my $if_report;

				if ($pileup_freq >= $report_cutoff || $tvc_freq >= $report_cutoff){
					$if_report = "Report";
				}else{
					$if_report = "NA";
				}

				my $if_PASS;
				if ($pileup_freq >= $report_cutoff and $tvc_freq >= $report_cutoff){
					$if_PASS = "PASS";
				}elsif ($pileup_freq < $report_cutoff and $tvc_freq < $report_cutoff){
					$if_PASS = "PASS";
				}else{
					$if_PASS = "Fail";
				}

				print "K03455.1\t$pos\t$ref\t$alt\t$info\t$pileup_freq\t$if_PASS\t$if_report\n";
			}
		}
	}
}
	
