use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw/$Bin/;
use List::Util qw(sum);


# 2022-8-3
# longfei.fu

my ($freebayes_var_file,$drug_aa_list_file,$name,$outdir);

GetOptions(
	"f:s"   => \$freebayes_var_file,       # Needed
	"a:s"   => \$drug_aa_list_file,        # Needed
	"n:s"   => \$name,                     # Needed
	"od:s"  => \$outdir,                   # Needed
	) or die "unknown args\n";


my $annot_outfile = "$outdir/$name\.freebayes.var.annot.xls";
open ANNOT, ">$annot_outfile" or die;
print ANNOT "Protein\tAA.Pos\tDNA.Pos\tRef.DNA\tRef.AA\tMut.AA.Annot\tMut.AA.Freq\n";


my %variants;
open VAR, "$freebayes_var_file" or die;
<VAR>;
while (<VAR>){
	chomp;
	my @arr = split /\t/;
	my $chr = $arr[0];
	my $pos = $arr[1];
	my $ref = $arr[2];
	my $alt = $arr[3];
	my $alt_depth = $arr[4];
	my $total_depth = $arr[5];
	my $alt_freq = $arr[6];

	my $var = "$chr\t$pos\t$ref\t$alt\t$alt_depth\t$total_depth\t$alt_freq";
	push @{$variants{$pos}}, $var;
}
close VAR;

# R41K
# Protein AA.Pos  DNA.Pos DNA     AA      AA.Long
# Protease        1       2253-2255       CCT     P       Proline
# Protease        2       2256-2258       CAG     Q       Glutamine
# Protease        3       2259-2261       GTC     V       Valine

open AA, "$drug_aa_list_file" or die;
<AA>;
while (<AA>){
	chomp;
	my @arr = split /\t/;
	my @region = split /\-/, $arr[2];
	
	my $first_pos = $region[0];
	my $second_pos = $first_pos + 1;
	my $third_pos = $region[1];
	
	my @pos;
	push @pos, $first_pos;
	push @pos, $second_pos;
	push @pos, $third_pos;

	my $ref_codon = $arr[-3];
	my @ref_base = split //, $arr[-3];

	my @first_pos_base; # include ref base
	my @second_pos_base;
	my @third_pos_base;

	push @first_pos_base, $ref_base[0];
	push @second_pos_base, $ref_base[1];
	push @third_pos_base, $ref_base[2];

	# Protease        2       2256-2258       CAG     Q       Glutamine
	# K03455.1        2258    G       A       6303    6304    1.00    GG:AA:2258

	# 2256: [C]
	# 2257: [A]
	# 2258: [G,A]



	# Protease        94      2532-2534       GGT     G       Glycine
	# K03455.1        2534    T       G       1583    3436    0.46
	# K03455.1        2534    T       C       738     3436    0.21

	# 2532: [G]
	# 2533: [G]
	# 2534: [T,G,C]

	# Protease        61      2433-2435       CAG     Q       Glutamine
	# K03455.1        2433    C       G       6782    7081    0.96    TCAGATACTC:AGAAATAGCT:2432
	# K03455.1        2435    G       A       6782    7081    0.96    TCAGATACTC:AGAAATAGCT:2432

	# 2433: [C,G]
	# 2434: [A]
	# 2435: [G,A]

	# CAG,CAA,GAG,GAA


	my %first_pos_alt_freq;
	my %second_pos_alt_freq;
	my %third_pos_alt_freq;

	my $idx = 0;
	for my $pos (@pos){
		$idx += 1;
		if (exists $variants{$pos}){
			my @var = @{$variants{$pos}};
			if ($idx == 1){
				for my $var (@var){
					# "$chr\t$pos\t$ref\t$alt\t$alt_depth\t$total_depth\t$alt_freq";
					my @val = split /\t/, $var;
					my $alt_base = $val[3];
					my $alt_base_freq = $val[-1];

					push @first_pos_base, $alt_base;
					$first_pos_alt_freq{$alt_base} = $alt_base_freq;
				}
			}elsif ($idx == 2){
				for my $var (@var){
					my @val = split /\t/, $var;
					my $alt_base = $val[3];
					my $alt_base_freq = $val[-1];

					push @second_pos_base, $alt_base;
					$second_pos_alt_freq{$alt_base} = $alt_base_freq;
			}else{
				# idx == 3
				for my $var (@var){
					my @val = split /\t/, $var;
					my $alt_base = $val[3];
					my $alt_base_freq = $val[-1];

					push @third_pos_base, $alt_base;
					$third_pos_alt_freq{$alt_base} = $alt_base_freq;
			}
		}
	}


	for my $b1 (@first_pos_base){
		for my $b2 (@second_pos_base){
			for my $b3 (@third_pos_base){
				my @nnn;
				push @nnn, $b1;
				push @nnn, $b2;
				push @nnn, $b3;

				my $codon = join("",@nnn);

				next if ($codon eq $ref_codon);

				my $mut_aa = ;

				my @aa_freq;

				my $b1_freq;
				if ($b1 ne $ref_base[0]){
					# alt base
					$b1_freq = $first_pos_alt_freq{$b1};
				}else{
					$b1_freq = 0;
				}

				my $b2_freq;
				if ($b2 ne $ref_base[1]){
					# alt base
					$b2_freq = $second_pos_alt_freq{$b2};
				}else{
					$b2_freq = 0;
				}

				my $b3_freq;
				if ($b3 ne $ref_base[2]){
					# alt base
					$b3_freq = $third_pos_alt_freq{$b3};
				}else{
					$b3_freq = 0;
				}

				push @aa_freq, $b1_freq;
				push @aa_freq, $b2_freq;
				push @aa_freq, $b3_freq;

				my @eff_freq;
				for my $freq (@aa_freq){
					if ($freq != 0){
						push @eff_freq, $freq;
					}
				}

				my $freq_sum;
				for my $freq (@eff_freq){
					$freq_sum += $freq;
				}

				my $alt_pos_num = scalar(@eff_freq);

				my $aa_freq = sprintf "%.2f", $freq_sum / $alt_pos_num;

				my $aa_var = $arr[-2].$arr[1].$mut_aa; # R41K
				print ANNOT "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t$aa_var\t$aa_freq\n";
			}
		}
	}
}
close AA;
close ANNOT;


