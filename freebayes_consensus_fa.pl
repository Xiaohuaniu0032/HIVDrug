use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw/$Bin/;
use List::Util qw(sum);
use POSIX qw(ceil);
use POSIX qw(floor);


my ($gvcf,my $ref_fa,$target_bed,$depth_cutoff,$freq_cutoff,$outdir);

GetOptions(
	"gvcf:s"    => \$gvcf,           # NEED
	"ref:s"     => \$ref_fa,         # NEED
	"bed:s"     => \$target_bed,     # NEED
	"d:i"       => \$depth_cutoff,   # Default: 20
	"f:f"       => \$freq_cutoff,    # Default: 0.05 (5%)
	"outdir:s"  => \$outdir,         # NEED
    ) or die "unknown args\n";


if (not defined $depth_cutoff){
	$depth_cutoff = 20;
}

if (not defined $freq_cutoff){
	$freq_cutoff = 0.05;
}


print "############ Analysis Param ############\n";
print "gvcf file is: $gvcf\n";
print "bed file is: $target_bed\n";
print "depth_cutoff is: $depth_cutoff [default: 20]\n";
print "freq_cutoff is: $freq_cutoff [default: 0.05 (5\%\)]\n";
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

sub get_ref_base{
	my ($target,$ref) = @_;

}

my $del_cutoff = 50; # 50%
my $ins_cutoff = 50; # 50%

# the first variant's pos

# the last variants's pos

my @cons_fa;

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
	
	if ($alt_allele =~ /\*/){
		# ref allele
		my @info = split /\;/, $arr[-3]; # DP=5045;END=2041;MIN_DP=5017
		my $min_depth = (split /\=/, $info[2])[1]; # 5017
		my $end_pos = (split /\=/, $info[1])[1]; # 2041
		my $target = "$ref_name\:$pos\-$end_pos";
		my $ref_base = &get_ref_base($target,$ref_fa);
		chomp $ref_base;
		push @cons_fa, $ref_base;
	}else{
		my @gt_info = split /\:/, $arr[-1]; # GT:DP:AD:RO:QR:AO:QA:GL => 0/0:5005:4752,171:4752:104017:171:1740:0,-1326.07,-8172.33
		my $gt = $gt_info[0];
		my $depth = $gt_info[1];
		my @allele_depth = split /\,/, $gt_info[2]; # 4752,171
		if ($alt_allele =~ /\,/){
			# multi alt allele
			my @alt_allele = split /\,/, $alt_allele;
			for my $alt_allele (@alt_allele){
				my $alt_allele_len = length($alt_allele);
				if ($ref_allele_len == $alt_allele_len){
					# SNP/MNP
					
				}

			}
		}else{
			# single alt allele
		}

	}
}
close GVCF;




