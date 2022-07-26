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

my ($gvcf,$name,$ref_fa,$target_bed,$samtools_bin,$depth_cutoff,$snp_freq_cutoff,$indel_freq_cutoff,$outdir);

GetOptions(
	"gvcf:s"        => \$gvcf,               # NEED
	"n:s"           => \$name,               # NEED
	"ref:s"         => \$ref_fa,             # NEED
	"bed:s"         => \$target_bed,         # NEED
	"samtools:s"    => \$samtools_bin,       # Default: /usr/bin/samtools
	"d:i"           => \$depth_cutoff,       # Default: 20
	"snp_f:f"       => \$snp_freq_cutoff,    # Default: 0.05 (5%)
	"indel_f:f"     => \$indel_freq_cutoff,  # Default: 0.6 (60%)
	"outdir:s"      => \$outdir,             # NEED
    ) or die "unknown args\n";


# check needed args
if (not defined $gvcf || not defined $name || not defined $ref_fa || not defined $target_bed || not defined $outdir){
	die "please check your args. [-gvcf, -n, -ref, -bed, -outdir are needed!]\n";
}

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

# check freq limit
if ($snp_freq_cutoff <= 0 || $snp_freq_cutoff >= 1 ){
	die "-snp_f limit is (0,1). default is: 0.05 (5\%)\n";
}

if ($indel_freq_cutoff <= 0 || $indel_freq_cutoff >= 1){
	die "-indel_f limit is (0,1). default is: 0.6 (60\%)\n";
}

my $freq_int = $snp_freq_cutoff * 100; # 0.05 => 5
my $freq_new = "Freq".$freq_int; # Freq5
my $consensus_fasta = "$outdir/$name\.freebayesConsensus\.$freq_new\.fasta";
my $log = "$outdir/$name\.$freq_new\.log";
my $var_summary = "$outdir/$name\.freebayes.variants\.$freq_new\.summary";


open LOG, ">$log" or die;
open SUMMARY, ">$var_summary" or die;

open CONS, ">$consensus_fasta" or die;
my $header = "$name\.freebayesConsensus\.$freq_new"; # IonCode_0101.freebayesConsensus.Freq20
print CONS "\>$header\n";

print LOG "############ Analysis Param ############\n";
print LOG "gvcf file is: $gvcf\n";
print LOG "bed file is: $target_bed\n";
print LOG "depth_cutoff is: $depth_cutoff [default: 20]\n";
print LOG "snp_freq_cutoff is: $snp_freq_cutoff [default: 0.05]\n";
print LOG "indel_freq_cutoff is: $indel_freq_cutoff [default: 0.6]\n";
print LOG "output dir is: $outdir\n"; 
print LOG "consensus.fasta file is: $consensus_fasta\n";
print LOG "log file is: $log\n";
print LOG "Variants summary file is: $var_summary\n";

my $ref_name = "K03455.1";


my %iupac_code_dict; # https://www.bioinformatics.org/sms/iupac.html
$iupac_code_dict{"A"} = "A";
$iupac_code_dict{"T"} = "T";
$iupac_code_dict{"C"} = "C";
$iupac_code_dict{"G"} = "G";

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

$iupac_code_dict{"AAA"} = "A";
$iupac_code_dict{"TTT"} = "T";
$iupac_code_dict{"CCC"} = "C";
$iupac_code_dict{"GGG"} = "G";

$iupac_code_dict{"CGT"} = "B";
$iupac_code_dict{"CTG"} = "B";
$iupac_code_dict{"GCT"} = "B";
$iupac_code_dict{"GTC"} = "B";
$iupac_code_dict{"TCG"} = "B";
$iupac_code_dict{"TGC"} = "B";


$iupac_code_dict{"AGT"} = "D";
$iupac_code_dict{"ATG"} = "D";
$iupac_code_dict{"GAT"} = "D";
$iupac_code_dict{"GTA"} = "D";
$iupac_code_dict{"TAG"} = "D";
$iupac_code_dict{"TGA"} = "D";


$iupac_code_dict{"ACT"} = "H";
$iupac_code_dict{"ATC"} = "H";
$iupac_code_dict{"CAT"} = "H";
$iupac_code_dict{"CTA"} = "H";
$iupac_code_dict{"TAC"} = "H";
$iupac_code_dict{"TCA"} = "H";


$iupac_code_dict{"ACG"} = "V";
$iupac_code_dict{"AGC"} = "V";
$iupac_code_dict{"CAG"} = "V";
$iupac_code_dict{"CGA"} = "V";
$iupac_code_dict{"GAC"} = "V";
$iupac_code_dict{"GCA"} = "V";


my @cons_fa;
my @processed_var;

my @low_depth_var;
my @snp_keep;
my @snp_skip;
my @indel_keep;
my @indel_skip;

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
	my $processed_var_num = scalar(@processed_var);

	#print "$_\n";
	#print "input_cons_fa: @cons_fa\n";
	if ($alt_allele =~ /\*/){
		# ref allele
		if ($processed_var_num == 0){
			# 第一行
			&process_ref_line($_,\@cons_fa,\@processed_var);
		}else{
			# 检查前一个变异位点
			my $former_var_end_pos = &get_former_end_pos(\@processed_var);
			my $exp_pos = $former_var_end_pos + 1;
			if ($pos == $exp_pos){
				# 位置正确
				&process_ref_line($_,\@cons_fa,\@processed_var);
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
					print LOG "gap: $target\n";
					my $ref_base = &get_ref_base($target,$ref_fa);
					push @cons_fa, $ref_base;
				}else{
					# 如果前一个变异最右端位置cover到了部分或者全部当前的变异,该如何处理
					# 暂时没有发现这种情况
				}
				&process_ref_line($_,\@cons_fa,\@processed_var); # 处理完gap碱基后,处理当前VCF行
			}
		}
	}else{
		# alt allele
		my @pass_allele; # v2版本: 如果ref allele的频率满足阈值,也考虑进来

		if ($processed_var_num == 0){
			# 第一行
			&process_var_line($_,\@pass_allele,\@processed_var,\@cons_fa);
		}else{
			# 检查前一个变异位点
			my $former_var_end_pos = &get_former_end_pos(\@processed_var);
			my $exp_pos = $former_var_end_pos + 1;
			if ($pos == $exp_pos){
				# 位置正确
				&process_var_line($_,\@pass_allele,\@processed_var,\@cons_fa);
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
				&process_var_line($_,\@pass_allele,\@processed_var,\@cons_fa);
			}
		}
	}
}
			
close GVCF;
for my $base (@cons_fa){
	print LOG "$base\n";
}
#print LOG "@cons_fa\n";
my $cons_fa = join("", @cons_fa);
print CONS "$cons_fa\n";
close CONS;
close LOG;



print SUMMARY "Low Depth Variants:\n";
for my $var (@low_depth_var){
	print SUMMARY "\t$var\n";
}
print SUMMARY "\n\n";

print SUMMARY "SNP Keep Variants:\n";
for my $var (@snp_keep){
	print SUMMARY "\t$var\n";
}
print SUMMARY "\n\n";

print SUMMARY "INDEL Keep Variants:\n";
for my $var (@indel_keep){
	print SUMMARY "\t$var\n";
}
print SUMMARY "\n\n";

print SUMMARY "SNP Skip Variants:\n";
for my $var (@snp_skip){
	print SUMMARY "\t$var\n";
}
print SUMMARY "\n\n";

print SUMMARY "INDEL Skip Variants:\n";
for my $var (@indel_skip){
	print SUMMARY "\t$var\n";
}
print SUMMARY "\n\n";




##################################################################
########################## sub function ##########################
##################################################################

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


sub get_former_end_pos{
	my ($processed_var_aref) = @_;
	my @var_list = @{$processed_var_aref};
	my $last_var = pop @var_list;
	#print "last_var: $last_var\n";
	my @var = split /\:/, $last_var; # chr/pos/ref/alt/...
	my $start_pos = $var[1];
	my $ref_allele = $var[2];
	my $ref_allele_len = length($ref_allele);
	my $end_pos = $start_pos + $ref_allele_len - 1;

	return($end_pos);
}

sub process_ref_line{
	my ($vcf_line,$cons_fa_aref,$processed_var_aref) = @_;
	my @v = @{$cons_fa_aref};
	my @arr = split /\t/, $vcf_line;
	my $pos = $arr[1]; # 2024
	#my $ref_allele = $arr[3];
	my @info = split /\;/, $arr[-3]; # DP=5045;END=2041;MIN_DP=5017
	my $depth = (split /\=/, $info[0])[1]; # 不适合用MIN_DP进行过滤
	my $end_pos = (split /\=/, $info[1])[1]; # 2041
	my $target = "$ref_name\:$pos\-$end_pos";
	my $ref_allele = &get_ref_base($target,$ref_fa);

	if ($depth >= $depth_cutoff){
		my $ref_base = &get_ref_base($target,$ref_fa);
		chomp $ref_base;
		push @{$cons_fa_aref}, $ref_base;
	}else{
		my $N_base = 'N' x ($end_pos - $pos + 1);
		push @{$cons_fa_aref}, $N_base;
	}

	push @{$processed_var_aref}, "$ref_name\:$pos\:$ref_allele";
}

sub process_var_line{
	my ($vcf_line,$pass_allele_aref,$processed_var_aref,$cons_fa_aref) = @_;
	#print "$vcf_line\n";
	my @arr = split /\t/, $vcf_line;
	my $pos = $arr[1];

	my @gt_info = split /\:/, $arr[-1]; # GT:DP:AD:RO:QR:AO:QA:GL => 0/0:5005:4752,171:4752:104017:171:1740:0,-1326.07,-8172.33
	
	my $total_depth = $gt_info[1];
	
	my $ref_allele = $arr[3];
	my $alt_allele = $arr[4];
	
	my @allele_depth = split /\,/, $gt_info[2]; # 4752,171 [ref/alt-1/alt-2/...]
	
	my $ref_allele_len = length($ref_allele);

	if ($total_depth < $depth_cutoff){
		my $N_base = 'N' x ($ref_allele_len);
		push @{$cons_fa_aref}, $N_base;
		my $var = "$ref_name\:$pos\:$ref_allele\:$alt_allele\:Depth\:$total_depth";
		push @{$processed_var_aref}, $var;
		print LOG "[Skip this Variant: Low Depth] $var\n";
		push @low_depth_var, $var;
	}else{
		my @allele;
		push @allele, $ref_allele;
		
		if ($alt_allele =~ /\,/){
			my @val = split /\,/, $alt_allele;
			for my $v (@val){
				push @allele, $v;
			}
		}else{
			push @allele, $alt_allele;
		}


		my %pass_allele_depth; # 记录pass allele深度信息

		my $idx = 0;
		for my $allele (@allele){
			if ($allele eq $ref_allele){
				# 检查ref allele频率
				my $ref_allele_depth = $allele_depth[0];
				my $ref_allele_freq;
				if ($ref_allele_depth > 0 and $total_depth > 0){
					$ref_allele_freq = sprintf "%.3f", $ref_allele_depth / $total_depth;
				}else{
					$ref_allele_freq = 0;
				}

				if ($ref_allele_freq >= $snp_freq_cutoff){
					push @{$pass_allele_aref}, $ref_allele;
					$pass_allele_depth{$ref_allele} = $ref_allele_depth;
				}
			}else{
				$idx += 1;
				my $alt_allele_depth = $allele_depth[$idx];
				my $alt_allele_len = length($allele);

				my $alt_allele_freq;
				if ($alt_allele_depth > 0 and $total_depth > 0){
					$alt_allele_freq = sprintf "%.3f", $alt_allele_depth / $total_depth;
				}else{
					$alt_allele_freq = 0;
				}

				my $var = "$ref_name\:$pos\:$ref_allele\:$allele\:$alt_allele_freq";

				if ($ref_allele_len == $alt_allele_len){
					# SNP/MNP
					# check freq
					if ($alt_allele_freq >= $snp_freq_cutoff){
						print LOG "[PASS, Written in Consensus Fasta] $var\n";
						push @{$pass_allele_aref}, $allele;
						$pass_allele_depth{$allele} = $alt_allele_depth;
						my $VAR = "$ref_name\:$pos\:$ref_allele\:$allele\:Freq\=$alt_allele_freq";
						push @snp_keep, $VAR;
					}else{
						# 频率不满足
						print LOG "[Skip this Variant: Low Freq] $var\n";
						my $VAR = "$ref_name\:$pos\:$ref_allele\:$allele\:Freq\=$alt_allele_freq";
						push @snp_skip, $VAR;
					}
				}else{
					# INDEL
					# 检查频率
					if ($alt_allele_freq >= $indel_freq_cutoff){
						print LOG "[Warning: INDEL!] [PASS, Written in Consensus Fasta] $var\n";
						#push @{$pass_allele_aref}, $allele;
						#$pass_allele_depth{$allele} = $alt_allele_depth;
						my $VAR = "$ref_name\:$pos\:$ref_allele\:$allele\:Freq\=$alt_allele_freq";
						push @indel_keep, $VAR;
					}else{
						# 频率不满足
						print LOG "[Warning: INDEL!] [Skip this Variant: Low Freq] $var\n";
						my $VAR = "$ref_name\:$pos\:$ref_allele\:$allele\:Freq\=$alt_allele_freq";
						push @indel_skip, $VAR;
					}
				}
			}
		}
		push @{$processed_var_aref}, "$ref_name\:$pos\:$ref_allele\:$alt_allele";

		
		# 处理pass_allele,判断哪些passed allele可以写入cons fa
		my $pass_allele_num = scalar(@{$pass_allele_aref});

		if ($pass_allele_num == 0){
			# no ref / no alt allele pass freq cutoff
			# may not happen
			my $end_pos = $pos + length($ref_allele) - 1;
			my $target = "$ref_name\:$pos\-$end_pos";
			my $ref_base = &get_ref_base($target,$ref_fa);
			push @{$cons_fa_aref}, $ref_base;
		}elsif ($pass_allele_num == 1){
			push @{$cons_fa_aref}, $pass_allele_aref->[0];
		}elsif ($pass_allele_num == 2){
			# 有2个满足条件的allele
			# 两个allele是否都是MNP
			my $allele_1_len = length($pass_allele_aref->[0]);
			my $allele_2_len = length($pass_allele_aref->[1]);
			
			if ($allele_1_len == $allele_2_len){
				# MNP
				my $idx = 0;
				for my $pos (1..$allele_1_len){
					#my @base;
					my %base;
					for my $allele (@{$pass_allele_aref}){
						my @allele_base = split //, $allele;
						my $base = $allele_base[$idx];
						#push @base,$base;
						$base{$base} = 1;
					}
					$idx += 1;

					my @base_str = keys %base;
					my $key = join("", @base_str);
					my $BASE = $iupac_code_dict{$key};
					push @{$cons_fa_aref}, $BASE;
				}
			}else{
				# 判断哪个allele的reads多,选择哪个allele生成一致性序列
				# 和generateConsensus一致
				my $first_allele = $pass_allele_aref->[0];
				my $second_allele = $pass_allele_aref->[1];

				my $first_allele_depth = $pass_allele_depth{$first_allele};
				my $second_allele_depth = $pass_allele_depth{$second_allele};

				if ($first_allele_depth >= $second_allele_depth){
					push @{$cons_fa_aref}, $first_allele;
				}else{
					push @{$cons_fa_aref}, $second_allele;
				}
			}
		}else{
			my @allele_sort_by_cov;
			foreach my $allele (sort {$pass_allele_depth{$b} <=> $pass_allele_depth{$a}} keys %pass_allele_depth){
				push @allele_sort_by_cov, $allele;
			}
			
			my @top_3_allele;
			push @top_3_allele, $allele_sort_by_cov[0];
			push @top_3_allele, $allele_sort_by_cov[1];
			push @top_3_allele, $allele_sort_by_cov[2];

			#print "pass_allele_is: @top_3_allele\n";
			
			my %PASS_ALLELE;
			for my $allele (@top_3_allele){
				$PASS_ALLELE{$allele} = 1;
			}
			
			my $check_alt_allele_len = 0;
			for my $allele (@top_3_allele){
				next if ($allele eq $ref_allele);
				my $len = length($allele);
				if ($len =! $ref_allele_len){
					$check_alt_allele_len = 1;
					last;
				}
			}
			#print "check_alt_allele_len: $check_alt_allele_len\n";

			if ($check_alt_allele_len == 0){
				# 所有alt allele长度与ref allele长度相同
				my $idx = 0;
				for my $pos (1..$ref_allele_len){
					#my @base;
					my %base;
					for my $allele (@top_3_allele){
						my @allele_base = split //, $allele;
						my $base = $allele_base[$idx];
						#push @base,$base;
						#print "base_is:$base\n";
						$base{$base} = 1;
					}
					$idx += 1;

					my @base_str = keys %base;
					my $key = join("",@base_str);
					my $BASE = $iupac_code_dict{$key};
					push @{$cons_fa_aref}, $BASE;
				}
			}else{
				# 哪个allele cov高选哪个allele
				# 和generateConsensus保持一致
				my @allele_sort_by_cov;

				foreach my $allele (sort {$pass_allele_depth{$b} <=> $pass_allele_depth{$a}} keys %pass_allele_depth){
					push @allele_sort_by_cov, $allele;
				}	

				my $top1_cov_allele = $allele_sort_by_cov[0];

				push @{$cons_fa_aref}, $top1_cov_allele;
			}
		}
	}
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
