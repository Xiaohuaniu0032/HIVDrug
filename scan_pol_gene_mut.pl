use strict;
use warnings;
use Data::Dumper;
use FindBin qw/$Bin/;

my ($pileup_var_file,$outfile) = @ARGV;

my $freq_cutoff = 0.05; # 5% mut freq cutoff
my $depth_cutoff = 100; # total depth cutoff
my $drug_aa_list = "/data/fulongfei/git_repo/HIVDrug/drug.aa.list.xls";



my %var;
open IN, "$pileup_var_file" or die;
my $header = <IN>;
chomp $header;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $pos = $arr[2];
	my $ref_base = $arr[3];
	my $alt_base = $arr[4];
	my $depth = $arr[-2];
	my $freq = $arr[-1];
	
	next if ($depth < $depth_cutoff);
	next if ($freq < $freq_cutoff);
	next if (length($ref_base) != length($alt_base)); # only keep snp

	push @{$var{$pos}}, "$ref_base\t$alt_base\t$freq";
}
close IN;


my %codon_list_short;
$codon_list{"ATT"} = "I"; $codon_list{"ATC"} = "I"; $codon_list{"ATA"} = "I";
$codon_list{"CTT"} = "L"; $codon_list{"CTC"} = "L"; $codon_list{"CTA"} = "L"; $codon_list{"CTG"} = "L"; $codon_list{"TTA"} = "L"; $codon_list{"TTG"} = "L";
$codon_list{"GTT"} = "V"; $codon_list{"GTC"} = "V"; $codon_list{"GTA"} = "V"; $codon_list{"GTG"} = "V";
$codon_list{"TTT"} = "F"; $codon_list{"TTC"} = "F";
$codon_list{"ATG"} = "M";
$codon_list{"TGT"} = "C"; $codon_list{"TGC"} = "C";
$codon_list{"GCT"} = "A"; $codon_list{"GCC"} = "A"; $codon_list{"GCA"} = "A"; $codon_list{"GCG"} = "A";




open AA, "$drug_aa_list" or die;
<AA>;
while (<AA>){
	chomp;
	my @arr = split /\t/;
	my @pos = split /\-/, arr[2];
	my $first = $pos[0];
	my $second = $first + 1;
	my $third = $pos[1];

	my $nt_three = $arr[3]; # CCT
	my @nt_three = split //, $nt_three;

	undef @pos;
	push @pos, $first;
	push @pos, $second;
	push @pos, $third;

	my @aaa_nt;
	push @aaa_nt, $nt_three[0];
	my @bbb_nt;
	push @bbb_nt, $nt_three[1];
	my @ccc_nt;
	push @ccc_nt, $nt_three[2];

	my $idx = 0;
	for my $pos (@pos){
		$idx += 1;
		if (exists $var{$pos}){
			my @var = @{$var{$pos}};
			for my $var (@var){
				my $alt_base = (split /\t/, $var)[1];
				if ($idx == 1){
					push @aaa_nt, $alt_base;
				}
				if ($idx == 2){
					push @bbb_nt, $alt_base;
				}
				if ($idx == 3){
					push @ccc_nt, $alt_base;
				}
			}
		}
	}

	my @var_BASE;
	for my $base_A (@aaa_nt){
		for my $base_B (@bbb_nt){
			for my $base_C (@ccc_nt){
				my $BASE = $base_A.$base_B.$base_C;
				next if ($BASE eq $nt_three); # skip ref
				push @var_BASE, $BASE;
			}
		}
	}

	# IonCode_0101    K03455.1        2258    G       A       6303    6304    1.0
	# IonCode_0101    K03455.1        2259    G       A       6320    6321    1.0
	# IonCode_0101    K03455.1        2282    C       T       1872    1876    0.998
	# IonCode_0101    K03455.1        2285    C       T       1872    1888    0.992
	# IonCode_0101    K03455.1        2288    A       C       4285    4302    0.996
	# IonCode_0101    K03455.1        2534    T       C       837     4066    0.206
	# IonCode_0101    K03455.1        2534    T       G       1651    4066    0.406


	# Protease        2       2256-2258       CAG     Q       Glutamine

	# 2256 [C]
	# 2257 [A]
	# 2258 [G,A]

	# CAG,CAA
	# ref is CAG, then mut aa is [CAA]




