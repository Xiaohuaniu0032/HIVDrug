use strict;
use warnings;
use File::Basename;
use FindBin qw/$Bin/;

# inpput file
# 1) drug.aa.list.xls
# Integrase       51      4380-4382       CAT     H       Histidine
# Integrase       66      4425-4427       ACA     T       Threonine
# Integrase       92      4503-4505       GAA     E       Glutamic acid
# Integrase       95      4512-4514       CAG     Q       Glutamine
# Integrase       97      4518-4520       ACA     T       Threonine

# 2) db/ScoresPI_1655774816824.tsv
# Rule    Position        AA      BIC     CAB     DTG     EVG     RAL
# H51Y    51      Y       10      15      10      15      15
# T66A    66      A       0       0       0       60      15
# T66I    66      I       5       10      5       60      15
# T66K    66      K       15      20      15      60      60
# E92G    92      G       0       10      0       30      15
# E92Q    92      Q       10      15      10      60      30
# E92V    92      V       0       10      0       60      30
# Q95K    95      K       0       0       0       10      10
# T97A    97      A       0       0       0       10      10


# Temp Format:
# {TAT,TAC} => Y
# Integrase H51Y 4380-4382 CAT TAT,TAC
# Integrase T66A,T66I,T66K 4425-4427 ACA GCT,GCC,GCA,GCG;ATT,ATC,ATA;AAA,AAG

# output format
# Chr/Pos/ID/Ref/Alt
# K03455.1/4380/Integrase,H51Y/C/T
# K03455.1/4382/Integrase,H51Y/T/C
# K03455.1/4425/Integrase,T66A,T66I,T66K/A/G
# K03455.1/4426/Integrase,T66A,T66I,T66K/C/T,A
# K03455.1/4427/Integrase,T66A,T66I,T66K/A/T,C,G

my ($drug_aa_list_file,$score_tsv_dir,$outdir) = @ARGV;

# ScoresINSTI_1655774833308.tsv     # Integrase
# ScoresNNRTI_1655774829129.tsv     # RT
# ScoresNRTI_1655774824871.tsv      # RT
# ScoresPI_1655774816824.tsv        # Protease


my %drug_aa;
open IN, "$drug_aa_list_file" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr  = split /\t/;
	my $prot = $arr[0];  # Protease
	my $pos  = $arr[1];
	$drug_aa{$prot}{$pos} = "$arr[2]\t$arr[3]"; # Protease => {41 => 2670-2672 \t ATG}
}
close IN;


my %AA2DNA;
my $codon_file = "/data/fulongfei/git_repo/HIVDrug/codon.list";
open IN, "$codon_file" or die;
while (<IN>){
	chomp;
	next if /^\#/;
	next if /^Amino/; # skip header
	my @arr = split /\t/;
	my $aa = $arr[1];
	my $codon = $arr[-1];
	my @codon = split /\,/, $codon;
	my @val;
	for my $v (@codon){
		$v =~ s/^\s+//; # 删除开始的空白
		push @val, $v;
	}
	my $val = join(",", @val);
	$AA2DNA{$aa} = $val;
}
close IN;


my @tsv = glob "$score_tsv_dir/*.tsv";
my %aa_tsv;

for my $file (@tsv){
	#print "$file\n";
	my $name_tmp = basename($file);
	my $name     = (split /\_/, $name_tmp)[0]; # ScoresINSTI

	my $p;
	if ($name eq "ScoresINSTI"){
		$p = "Integrase";
	}elsif ($name eq "ScoresNNRTI"){
		$p = "RT";
	}elsif ($name eq "ScoresNRTI"){
		$p = "RT";
	}elsif ($name eq "ScoresPI"){
		$p = "Protease";
	}else{
		$p = "NA";
	}

	if ($p eq "NA"){
		die "please check your drug resistance tsv file\n";
	}

	open IN, "$file" or die;
	<IN>;
	while (<IN>){
		chomp;
		my @arr = split /\t+/;
		my $aa = $arr[0];
		next if ($aa =~ /\+/);
		next if ($aa =~ /del/ || $aa =~ /ins/); # only for SNP
		my $pos = $arr[1]; # 41
		push @{$aa_tsv{$p}{$pos}}, $aa; # RT => {65 => [K65E,K65N]}
	}
	close IN;
}

my $temp_file = "$outdir/HS.temp.txt";
open TEMP, ">$temp_file" or die;

my @prot = qw/Protease RT Integrase/;
for my $p (@prot){
	my @pos = sort ({$a <=> $b} keys %{$aa_tsv{$p}}); # POS:1-N
	for my $pos (@pos){
		# 每个位置可能有多个突变AA
		my $aa_var_aref = $aa_tsv{$p}{$pos};
		my @var = @{$aa_var_aref};
		
		my %var_dna_seq; # K65E, E : GAA,GAG
		# 得到突变氨基酸对应的DNA碱基
		for my $var (@var){
			# K65E
			# I need E
			my @val = split //, $var;
			my $var_aa = $val[-1]; # E
			my $codon = $AA2DNA{$var_aa}; # GAA,GAG
			my @codon = split /\,/, $codon;
			for my $c (@codon){
				$var_dna_seq{$c} = 1;
			}
		}

		my @var_dna_seq = keys %var_dna_seq;
		my $var_dna_seq = join(",",@var_dna_seq); # DNA碱基

		my $ref_info = $drug_aa{$p}{$pos};
		my $ref_pos = (split /\t/, $ref_info)[0]; # 参考基因组位置
		my $ref_seq = (split /\t/, $ref_info)[1]; # 参考基因组碱基
		my $val = join(",",@var); # 所有可能的唯一AA突变列表
		#print "$ref_info\n";
		print TEMP "$p\t$val\t$ref_pos\t$ref_seq\t$var_dna_seq\n";
	}
}
close TEMP;

my $hs_file = "$outdir/HIV.drug.HS.vcf";
open O, ">$hs_file" or die;
print O "\#CHROM\tPOS\tID\tREF\tALT\n";

# output format
# Chr/Pos/ID/Ref/Alt
# K03455.1/4380/Integrase,H51Y/C/T
# K03455.1/4382/Integrase,H51Y/T/C
# K03455.1/4425/Integrase,T66A,T66I,T66K/A/G
# K03455.1/4426/Integrase,T66A,T66I,T66K/C/T,A
# K03455.1/4427/Integrase,T66A,T66I,T66K/A/T,C,G


sub uniq_base{
	my ($ref_base,$alt_codon,$idx) = @_; # A; TTC,TTT
	my %uniq_base;
	if ($alt_codon =~ /\,/){
		my @alt_codon = split /\,/, $alt_codon;
		for my $v (@alt_codon){
			my @v = split //, $v;
			my $base = $v[$idx]; # 0-based
			$uniq_base{$base} += 1;
		}
	}
	my @base = keys %uniq_base;
	
	return(\@base);
}

open IN, "$temp_file" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $p = $arr[0];
	my @pos = split /\-/, $arr[2];
	my $start = $pos[0];
	my $end = $pos[1];
	
	my $ref = $arr[-2];
	my @ref = split //, $ref;

	my @alt = split /\,/, $arr[-1];

	my $ID = $arr[0].';'.$arr[1].';'.$arr[2].';'.$arr[3].';'.$arr[4];

	
	for my $idx (0..2){
		my $ref_base = $ref[$idx];
		my $ref_pos = $start + $idx;
		my $wait_base_aref = uniq_base($ref_base,$arr[-1],$idx);
		my @wait_base = @{$wait_base_aref};

		my @real_alt_base;
		for my $base (@wait_base){
			if ($base ne $ref_base){
				push @real_alt_base, $base;
			}
		}

		if (scalar(@real_alt_base) > 0){
			my $v = join(",", @real_alt_base);
			print O "K03455.1\t$ref_pos\t$ID\t$ref_base\t$v\n";
		}
	}
}
close IN;
close O;
