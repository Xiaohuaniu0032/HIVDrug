use strict;
use warnings;
use File::Basename;
use FindBin qw/$Bin/;


# input file: HIV.drug.HS.vcf

my ($infile) = @ARGV;
my $fa = "/data/fulongfei/git_repo/HIVDrug/ref/K03455.fasta";
my $samtools = "/usr/bin/samtools";
my $chr_name = "K03455.1";
print "track type=bedDetail name=\"HIVDurg_HotSpots\" description=\"HIV Drug Resistance Analysis\"\n";


#track type=bedDetail name="CHP2_HotSpots" description="CHP2_COSMIC_Mutations_v60"
#chr1    43814978        43814979        COSM27286       0       +       REF=G;OBS=A;ANCHOR=A    AMP195


open IN, "$infile" or die;
while (<IN>){
	chomp;
	next if (/^\#/);
	next if (/^$/);
	my @arr = split /\t/;
	my $pos = $arr[1]; # 1-based
	my $pos_former = $pos - 1;
	my $anchor_base_tmp = `$samtools faidx $fa $chr_name\:$pos_former\-$pos_former`; # >K03455.1:4924-4924\nA\n
	my $anchor_base     = (split /\n/, $anchor_base_tmp)[1];
	my $info = $arr[2]; # Protease;L10F;2280-2282;CTC;TTT,TTC
	my @info = split /\;/, $info;
	my $info_new = join("_",@info);

	my $alt = $arr[-1];
	if ($alt =~ /\,/){
		my @alt = split /\,/, $alt;
		for my $v (@alt){
			print "$chr_name\t$pos_former\t$pos\t$info_new\t0\t+\tREF\=$arr[-2]\;OBS\=$v\;ANCHOR\=$anchor_base\t$info_new\n";
		}
	}else{
		print "$chr_name\t$pos_former\t$pos\t$info_new\t0\t+\tREF\=$arr[-2]\;OBS\=$alt\;ANCHOR\=$anchor_base\t$info_new\n";
	}
}
close IN;