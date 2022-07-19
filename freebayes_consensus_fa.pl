use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw/$Bin/;
use List::Util qw(sum);

# 2022-7-20
# longfei.fu@thermofisher.com
# 601435543@qq.com

my ($gvcf,$ref_fa,$target_bed,$samtools_bin,$depth_cutoff,$snp_freq_cutoff,$indel_freq_cutoff,$outdir);

GetOptions(
	"gvcf:s"        => \$gvcf,               # NEED
	"ref:s"         => \$ref_fa,             # NEED
	"bed:s"         => \$target_bed,         # NEED
	"samtools:s"    => \$samtools_bin,       # Default: /usr/bin/samtools
	"d:i"           => \$depth_cutoff,       # Default: 20
	"snp_f:f"       => \$snp_freq_cutoff,    # Default: 0.05 (5%)
	"indel_f:f"     => \$indel_freq_cutoff,  # Default: 0.6 (60%)
	"outdir:s"      => \$outdir,             # NEED
    ) or die "unknown args\n";


if (not defined $depth_cutoff){
	$depth_cutoff = 20;
}

if (not defined $snp_freq_cutoff){
	$snp_freq_cutoff = 0.05; # 5%
}

if (not defined $indel_freq_cutoff){
	$indel_freq_cutoff = 0.6; # 60%
}

if (not defined $samtools_bin){
	$samtools_bin = "/usr/bin/samtools";
}

print "############ Analysis Param ############\n";
print "gvcf file is: $gvcf\n";
print "bed file is: $target_bed\n";
print "depth_cutoff is: $depth_cutoff [default: 20]\n";
print "snp_freq_cutoff is: $snp_freq_cutoff [default: 0.05]\n";
print "indel_freq_cutoff is: $indel_freq_cutoff [default: 0.6]\n";
print "output dir is: $outdir\n"; 


my $ref_name = "K03455.1";

# /data/fulongfei/analysis/hiv/JSCDC/IonCode_0101/freebayes/IonCode_0101.freebayes.ploidy2.gvcf
# /data/fulongfei/git_repo/HIVDrug/BED/new_BED/POL.bed
# Protease	1	2253-2255	CCT	P	Proline
# ///
# Protease	99	2547-2549	TTT	F	Phenylalanine
# RT	1	2550-2552	CCC	P	Proline
# ///
# RT	560	4227-4229	CTA	L	Leucine
# Integrase	1	4230-4232	TTT	F	Phenylalanine
# ///
# Integrase	288	5091-5093	GAT	D	Aspartic acid
# Integrase	289	5094-5096	TAG	Stop	Stop codons

my %iupac_code_dict; # https://www.bioinformatics.org/sms/iupac.html
$iupac_code_dict{"AA"} = "A";
$iupac_code_dict{"TT"} = "T";
$iupac_code_dict{"CC"} = "C";
$iupac_code_dict{"GG"} = "G";

$iupac_code_dict{"AG"} = "R";
$iupac_code_dict{"GA"} = "R";

$iupac_code_dict{"CT"} = "Y";
$iupac_code_dict{"TC"} = "Y";

$iupac_code_dict{"GC"} = "S";
$iupac_code_dict{"CG"} = "S";

$iupac_code_dict{"AT"} = "W";
$iupac_code_dict{"TA"} = "W";

$iupac_code_dict{"GT"} = "K";
$iupac_code_dict{"TG"} = "K";

$iupac_code_dict{"AC"} = "M";
$iupac_code_dict{"CA"} = "M";


my @cons_fa;

sub get_ref_base{
	my ($target,$ref) = @_;
	#fulongfei@ion-bfx:/data/fulongfei/git_repo/HIVDrug/ref$ samtools faidx K03455.fasta K03455.1:2024-2041
	#>K03455.1:2024-2041
	#GCTGTTGGAAATGTGGAA
	my $ref_base = (split /\n/, `$samtools_bin faidx $ref $target`)[1];
	chomp $ref_base;
	#print "ref base is: $ref_base\n";

	return $ref_base;
}



# Protease	1	2253-2255	CCT	P	Proline
# Protease	2	2256-2258	CAG	Q	Glutamine
# Protease	3	2259-2261	GTC	V	Valine
# ......
# Integrase	288	5091-5093	GAT	D	Aspartic acid
# Integrase	289	5094-5096	TAG	Stop	Stop codons

# K03455.1        2019    5259 => [2020,5259]

#my $start_pos = 2020;
#my $end_pos   = 5259;

# the first variant's pos
#my @first_line;
#open IN, "$gvcf" or die;
#while (<IN>){
#	chomp;
#	next if (/^\#/);
#	next if (/^$/);
#	if (/^K03455/){
#		push @first_line, $_;
#		last;
#	}
#}
#close IN;

#my $first_line = $first_line[0];
#undef @first_line;
#@first_line = split /\t/, $first_line;
#my $alt_allele = $first_line[4];
#if ($alt_allele =~ /\*/){
#	my $pos = $arr[1];
#	my $gap_len = abs($pos - $start_pos);
#	if ($gap_len >= 1){
#		if ($pos > $start_pos){
#			# add N at the head of consensus fasta
#
#		}
#	}
#}



# the last variants's pos
#my @last_line;
#my $total_line = (split /\s/, `wc -l $gvcf`)[0]; # 619 IonCode_0101.freebayes.ploidy2.gvcf
#my $line = 0;
#open IN, "$gvcf" or die;
#while (<IN>){
#	chomp;
#	$line += 1;
#	if ($line == $total_line){
#		push @last_line, $_;
#		last;
#	}
#}
#close IN;

sub get_former_end_pos{
	my ($processed_var_aref) = @_;
	my @var_list = @{$processed_var_aref};
	my $last_var = pop @var_list;
	my @var = split /\t/, $last_var; # pos/ref/alt
	my $start_pos = $var[0];
	my $ref_allele = $var[1];
	my $ref_allele_len = length($ref_allele);
	my $end_pos = $start_pos + $ref_allele_len - 1;

	return($end_pos);
}

sub process_ref_line{
	my ($vcf_line) = @_;
	my @arr = split /\t/, $vcf_line;
	my $pos = $arr[1];
	my @info = split /\;/, $arr[-3]; # DP=5045;END=2041;MIN_DP=5017
	my $depth = (split /\=/, $info[0])[1]; # 不适合用MIN_DP进行过滤
	my $end_pos = (split /\=/, $info[1])[1]; # 2041
	my $target = "$ref_name\:$pos\-$end_pos";
	
	if ($depth >= $depth_cutoff){
		my $ref_base = &get_ref_base($target,$ref_fa);
		chomp $ref_base;
		push @cons_fa, $ref_base;
	}else{
		my $N_base = 'N' x ($end_pos - $pos + 1);
		push @cons_fa, $N_base;
	}
}

sub process_var_line{
	my ($vcf_line) = @_;


}

sub process_pass_allele{
	my ($pass_allele_aref) = @_;

}

sub extend_base_for_left_cons{
	#
}

sub extend_base_for_right_cons{
	# 
}

# 4310位置后面应该紧跟4311,但VCF文件中下一个位置是从4312开始
# 这种特殊情况需要考虑,否则会出现遗漏 

#K03455.1        4305    .       T       <*>     0       .       DP=3911;END=4309;MIN_DP=3902    GQ:DP:MIN_DP:QR
#:RO:QA:AO        496066:3911:3902:498070:3897:2004:14
#K03455.1        4310    .       C       T       73012.9 .       AB=0;ABP=0;AC=2;AF=1;AN=2;AO=3739;CIGAR=1X;DP=3888;DPB=3888;DPRA=0;EPP=7.61052;EPPR=8.34296;GTI=0;LEN=1;MEANALT=3;MQM=46.3188;MQMR=19.3605;NS=1;NUMALT=1;ODDS=4361.16;PAIRED=0;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=95865;QR=4488;RO=147;RPL=3605;RPP=6999.95;RPPR=17.2061;RPR=134;RUN=1;SAF=1717;SAP=57.0358;SAR=2022;SRF=141;SRP=272.229;SRR=6;TYPE=snp;technology.IONTORRENT=1       GT:DP:AD:RO:QR:AO:QA:GL 1/1:3888:147,3739:147:4488:3739:95865:-8148.61,-954.567,0
#K03455.1        4312    .       T       <*>     0       .       DP=13509;END=4317;MIN_DP=3900   GQ:DP:MIN_DP:QR:RO:QA:AO        1.86135e+06:13509:3900:1872180:13437:10829:71

my @processed_var;

open GVCF, "$gvcf" or die;
while (<GVCF>){
	chomp;
	next if (/^\#/);
	next if (/^$/);
	my @arr = split /\t/;
	my $ref = $arr[0];
	my $pos = $arr[1];
	my $ref_allele = $arr[3];
	my $alt_allele = $arr[4];
	my $QUAL = $arr[5];

	my $ref_allele_len = length($ref_allele);
	
	my $var = "$pos\t$ref_allele\t$alt_allele";
	my $processed_var_num = scalar(@processed_var);

	if ($alt_allele =~ /\*/){
		# ref allele
		if ($processed_var_num == 0){
			# 第一行
			&process_ref_line($_);
		}else{
			# 检查前一个变异位点
			my $former_var_end_pos = &get_former_end_pos(\@processed_var);
			my $exp_pos = $former_var_end_pos + 1;
			if ($pos == $exp_pos){
				# 位置正确
				&process_ref_line($_);
			}else{
				if ($pos - $former_var_end_pos > 1){
					# 正常情况应该是相差1
					# 跳空的位置补上ref base
					my @gap_pos;
					my $sp = $former_var_end_pos + 1;
					my $ep = $pos - 1;
					for my $pos ($sp..$ep){
						push @gap_pos, $pos;
					}

					my $target = "$ref_name\:$gap_pos[0]\-$gap_pos[-1]";
					my $ref_base = &get_ref_base($target,$ref_fa);
					push @cons_fa, $ref_base;
				}else{
					# 如果前一个变异最右端位置cover到了部分或者全部当前的变异,该如何处理
					# 暂时没有发现这种情况
				}
				&process_ref_line($_); # 处理完gap碱基后,处理当前VCF行
			}
		}
	}else{
		# alt allele
		if ($processed_var_num == 0){
			# 第一行
			&process_var_line($_);
		}else{
			# 检查前一个变异位点
			my $former_var_end_pos = &get_former_end_pos(\@processed_var);
			my $exp_pos = $former_var_end_pos + 1;
			if ($pos == $exp_pos){
				# 位置正确
				&process_var_line($_);
			}else{
				if ($pos - $former_var_end_pos > 1){
					# 正常情况应该是相差1
					# 跳空的位置补上ref base



					push @cons_fa, $ref_base;
				}else{
					# 如果前一个变异最右端位置cover到了部分或者全部当前的变异,该如何处理
					# 暂时没有发现这种情况
				}
				&process_var_line($_);
			}
		}
		my @gt_info = split /\:/, $arr[-1]; # GT:DP:AD:RO:QR:AO:QA:GL => 0/0:5005:4752,171:4752:104017:171:1740:0,-1326.07,-8172.33
		my $gt = $gt_info[0];
		my $depth = $gt_info[1];

		if ($depth < $depth_cutoff){
			my $N_base = 'N' x ($ref_allele_len);
			push @cons_fa, $N_base;
			push @processed_var, $var;
			next;
		}


		my @allele_depth = split /\,/, $gt_info[2]; # 4752,171 (one alt allele, first '4752' is the ref allele depth) or 27,187,4821 (two alt alleles)
		
		my @pass_allele;

		# 第一步:生成pass_allele
		# 第二部:根据pass_allele,生成cons fa
		
		my %alt_allele_depth; # 记录alt allele的深度信息

		if ($alt_allele =~ /\,/){
			# multi alt allele
			my @alt_allele = split /\,/, $alt_allele;
			my $idx = 0;
			for my $allele (@alt_allele){
				$idx += 1;
				my $alt_allele_depth = $allele_depth[$idx];
				$alt_allele_depth{$allele} = $alt_allele_depth; # 记录alt allele的深度信息
				
				my $alt_allele_len = length($allele);

				my $alt_allele_freq;
				if ($alt_allele_depth > 0 and $depth > 0){
					$alt_allele_freq = sprintf "%.3f", $alt_allele_depth / $depth;
				}else{
					$alt_allele_freq = 0;
				}

				my $var = "$ref_name\:$pos\:$ref_allele\:$allele\:$alt_allele_freq";

				if ($ref_allele_len == $alt_allele_len){
					# SNP/MNP
					# check freq
					if ($alt_allele_freq >= $snp_freq_cutoff){
						print "[PASS, Written in Consensus Fasta] $var\n";
						push @pass_allele, $allele;
					}else{
						# 频率不满足
						print "[Skip this Variant: Low Freq] $var\n";
					}
				}else{
					# INDEL
					# 检查频率
					if ($alt_allele_freq >= $indel_freq_cutoff){
						print "[Warning: INDEL!] [PASS, Written in Consensus Fasta] $var\n";
						push @pass_allele, $allele;
					}else{
						# 频率不满足
						print "[Warning: INDEL!] [Skip this Variant: Low Freq] $var\n";
					}
				}
			}
		}else{
			# single alt allele
			my $alt_allele_depth = $allele_depth[1];
			my $alt_allele_len = length($alt_allele);
			
			$alt_allele_depth{$alt_allele} = $alt_allele_depth;

			# 突变频率
			my $alt_allele_freq;
			if ($alt_allele_depth > 0 and $depth > 0){
				$alt_allele_freq = sprintf "%.3f", $alt_allele_depth / $depth;
			}else{
				$alt_allele_freq = 0;
			}

			my $var = "$ref_name\:$pos\:$ref_allele\:$alt_allele\:$alt_allele_freq";

			# 检查是SNP/MNP还是INDEL
			if ($ref_allele_len == $alt_allele_len){
				# SNP/MNP
				# 检查频率
				if ($alt_allele_freq >= $snp_freq_cutoff){
					push @pass_allele, $alt_allele;
					print "[PASS, Written in Consensus Fasta] $var\n";
				}else{
					print "[Skip this Variant: Low Freq] $var\n";
				}
			}else{
				# INDEL
				# 检查频率
				if ($alt_allele_freq >= $indel_freq_cutoff){
					print "[Warning: INDEL!] [PASS, Written in Consensus Fasta] $var\n";
					push @pass_allele, $alt_allele;
				}else{
					print "[Warning: INDEL!] [Skip this Variant: Low Freq] $var\n";
				}
			}
		}

		# 处理pass allele
		my $alt_num = scalar(@pass_allele);
		
		if ($alt_num == 0){
			# no alt allele
			my $end_pos = $pos + length($ref_allele) - 1;
			my $target = "$ref_name\:$pos\-$end_pos";
			my $ref_base = &get_ref_base($target,$ref_fa);
			push @cons_fa, $ref_base;
		}elsif ($alt_num == 1){
			# 有一个满足条件的alt allele
			push @cons_fa, $pass_allele[0];
		}elsif ($alt_num == 2){
			# 有2个满足条件的alt allele
			# 两个allele是否都是MNP
			my $allele_1_len = length($pass_allele[0]);
			my $allele_2_len = length($pass_allele[1]);
			if ($allele_1_len == $allele_2_len){
				# MNP
				my $idx = 0;
				for my $pos (1..$allele_1_len){
					my @base;
					for my $allele (@pass_allele){
						my @allele_base = split //, $allele;
						my $base = $allele_base[$idx];
						push @base,$base;
					}
					$idx += 1;

					my $key = join("", @base);
					my $BASE = $iupac_code_dict{$key};
					push @cons_fa, $BASE;
				}
			}else{
				# 两个alt allele. 一个SNP 一个INDEL
				# 判断哪个allele的reads多,选择哪个allele生成一致性序列
				# 这种情况应该很少见
				my $first_allele = $pass_allele[0];
				my $second_allele = $pass_allele[1];

				my $first_allele_depth = $alt_allele_depth{$first_allele};
				my $second_allele_depth = $alt_allele_depth{$second_allele};

				if ($first_allele_depth >= $second_allele_depth){
					push @cons_fa, $first_allele;
				}else{
					push @cons_fa, $second_allele;
				}
			}
		}else{
			# 有>=3个满足条件的alt allele
			# 如果>=3个alt allele都是SNP,选择top2生成一致性序列
			# 检查多个allele长度是否都相同
			my $check_alt_allele_len = 0;
			for my $allele (@pass_allele){
				my $len = length($allele);
				if ($len =! $ref_allele_len){
					$check_alt_allele_len = 1;
					last;
				}
			}

			if ($check_alt_allele_len == 0){
				# 所有alt allele长度与ref allele长度相同
				# 取top2 depth alt allele
				my @allele_sort_by_cov;
				foreach my $allele (sort {$alt_allele_depth{$b} <=> $alt_allele_depth{$a}} keys %alt_allele_depth){
					push @allele_sort_by_cov, $allele;
				}

				my @top2_cov_allele;
				my $top1_allele = shift @allele_sort_by_cov;
				my $top2_allele = shift @allele_sort_by_cov;
				push @top2_cov_allele, $top1_allele;
				push @top2_cov_allele, $top2_allele;

				my $top1_allele_len = length($top1_allele);
				my $idx = 0;
				for my $pos (1..$top1_allele_len){
					my @base;
					for my $allele (@top2_cov_allele){
						my @allele_base = split //, $allele;
						my $base = $allele_base[$idx];
						push @base,$base;
					}
					$idx += 1;

					my $key = join("",@base);
					my $BASE = $iupac_code_dict{$key};
					push @cons_fa, $BASE;
				}
			}else{
				# >=3个alt allele中存在一个INDEL,则直接取cov最高的作为alt allele
				my @allele_sort_by_cov;
				foreach my $allele (sort {$alt_allele_depth{$b} <=> $alt_allele_depth{$a}} keys %alt_allele_depth){
					push @allele_sort_by_cov, $allele;
				}

				my $top1_cov_allele = $allele_sort_by_cov[0];

				push @cons_fa, $top1_cov_allele;
			}
		}
	}
	push @processed_var, $var;
}
			
close GVCF;
print "@cons_fa\n";
my $cons_fa = join("", @cons_fa);
print "Consensus Fasta is: $cons_fa\n";