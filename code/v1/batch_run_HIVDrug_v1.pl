use strict;
use warnings;
use File::Basename;
use FindBin qw/$Bin/;

# /data/zsq/HIV/Fst_BAMs
my ($bam_dir, $outdir) = @ARGV;
my @bam_files = glob "$bam_dir/*.bam";

for my $bam (@bam_files){
	my $base_name = basename($bam); # IonCode_0101_R_2022_07_05_22_42_19_user_GSS5-0655-1-20220705_20220707.bam
	my @name = split /\_/, $base_name;
	my $name = "$name[0]\_$name[1]";
	print "$name\n";

	`perl /data/fulongfei/git_repo/HIVDrug/code/v1/HIVDrug_v1.pl $bam $name $outdir`;
}

