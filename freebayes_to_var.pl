use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw/$Bin/;


my ($vcf_file,$name,$snp_cutoff,$indel_cutoff,$outdir);
my $ref_name = "K03455.1";

GetOptions(
	"vcf:s"     => \$vcf_file,          # Needed
	"n:s"       => \$name,              # Needed
 	"snp_f:f"   => \$snp_cutoff,        # Default: 0.05 (5%)
	"indel_f:f" => \$indel_cutoff,      # Default: 1 (100%)
	"od:s"      => \$outdir,            # Needed
	) or die "unknown args\n";


if (not defined $vcf_file || not defined $name || not defined $outdir){
	die "please check -vcf <vcf_file> -n <name> -od <outdir> args\n";
}

if (not defined $snp_cutoff){
	$snp_cutoff = 0.05;
}

if (not defined $indel_cutoff){
	$indel_cutoff = 1; # 100%
}

my $outfile = "$outdir/$name\.freebayes.var.freq.txt";
open O, ">$outfile" or die;
print O "CHROM\tPOS\tRef_Allele\tAlt_Allele\tAlt_Allele_Depth\tTotal_Depth\tAlt_Allele_Freq\tRef_Allele\:Pass_Allele\:Pos\n";

open VCF, "$vcf_file" or die;
while (<VCF>){
	chomp;
	next if (/^$/);
	next if (/^\#/);
	my @arr = split /\t/;
	my $chr = $arr[0];
	my $pos = $arr[1];
	my $ref_allele = $arr[3];
	my $alt_allele = $arr[4];
	my @info = split /\:/, $arr[-1]; # GT:DP:AD:RO:QR:AO:QA:GL => [0/0:5005:4752,171:4752:104017:171:1740:0,-1326.07,-8172.33]
	my $total_depth = $info[1];
	my @alt_allele_depth = split /\,/, $info[2];
	shift @alt_allele_depth; # remove ref allele depth

	my %alt_allele_depth;
	my @pass_allele;

	my $ref_allele_len = length($ref_allele);

	next if ($alt_allele =~ /\*/);

	if ($alt_allele =~ /\,/){
		my @alt_allele = split /\,/, $alt_allele;
		my $idx = 0;
		for my $alt_allele_tmp (@alt_allele){
			my $alt_allele_depth = $alt_allele_depth[$idx];
			$alt_allele_depth{$alt_allele_tmp} = $alt_allele_depth;

			my $alt_allele_len = length($alt_allele_tmp);

			if ($alt_allele_len == $ref_allele_len){
				my $freq;
				if ($alt_allele_depth > 0 and $total_depth > 0){
					$freq = sprintf "%.2f", $alt_allele_depth / $total_depth;
				}else{
					$freq = 0;
				}

				if ($freq >= $snp_cutoff){
					push @pass_allele, $alt_allele_tmp;
				}
				#print "$ref_allele\t$alt_allele_tmp\t$alt_allele_depth\t$total_depth\t$freq\n";
			}
			$idx += 1;
		}
	}else{
		my $alt_allele_len = length($alt_allele);

		if ($alt_allele_len == $ref_allele_len){
			my $alt_allele_depth = $alt_allele_depth[0];
			$alt_allele_depth{$alt_allele} = $alt_allele_depth;

			my $freq;
			if ($alt_allele_depth > 0 and $total_depth > 0){
				$freq = sprintf "%.2f", $alt_allele_depth / $total_depth;
			}else{
				$freq = 0;
			}

			if ($freq >= $snp_cutoff){
				push @pass_allele, $alt_allele;
			}
		}
	}

	my $pass_allele_num = scalar(@pass_allele);
	next if ($pass_allele_num == 0);
	#print "pass_allele_is: @pass_allele\n";

	my %alt_allele_depth_sum;
	# K03455.1        2110    .       TA      CT,CA
	# K03455.1        2110    T       C       547     6356    0.09    TA:CT
	# K03455.1        2110    T       C       5795    6356    0.91    TA:CA

	for my $alt_allele (@pass_allele){
		my @ref_base = split //, $ref_allele;
		my @alt_base = split //, $alt_allele;

		my $alt_allele_depth = $alt_allele_depth{$alt_allele};
		my $alt_allele_freq = sprintf "%.2f", $alt_allele_depth / $total_depth;

		my $note = "$ref_allele\:$alt_allele";
		
		my $base_Num = scalar(@ref_base);
		for my $i (1..$base_Num){
			my $idx = $i - 1;
			my $ref = $ref_base[$idx];
			my $alt = $alt_base[$idx];
			if ($ref ne $alt){
				my $POS = $pos + $idx;
				my $var = "$ref\t$alt"; # ref/alt
				$alt_allele_depth_sum{$POS}{$var} += $alt_allele_depth;
				#print O "$ref_name\t$POS\t$ref\t$alt\t$alt_allele_depth\t$total_depth\t$alt_allele_freq\t$note\n\n";
			}
		}
	}

	my $pass_allele = join(",",@pass_allele);
	my $VAR = "$ref_allele\:$pass_allele\:$pos";

	foreach my $pos (sort {$a <=> $b} keys %alt_allele_depth_sum){
		my @var = keys %{$alt_allele_depth_sum{$pos}};
		for my $var (@var){
			my $depth = $alt_allele_depth_sum{$pos}{$var};
			my $freq = sprintf "%.2f", $depth/$total_depth;
			print O "$ref_name\t$pos\t$var\t$depth\t$total_depth\t$freq\t$VAR\n";
		}
	}
}
close VCF;
close O;





