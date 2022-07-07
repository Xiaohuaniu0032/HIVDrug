use strict;
use warnings;
use File::Basename;

my ($outdir) = @ARGV;

my $resdir = "/data/fulongfei/analysis/hiv/generateCons_test";
my $tvc_freq_dir = "/data/fulongfei/analysis/hiv/re_analysis_new_ref";

my @parsed_files = glob "$resdir/*/*.parsed.txt";

for my $file (@parsed_files){
	my $name = (split /\./, basename($file))[0];
	my $tvc_freq_file = "$tvc_freq_dir/$name\_rawlib/freq_cmp.txt";
	my $outfile = "$outdir/$name/$name\.Cons.Check.txt";
	#print "$outfile\n";
	`perl /data/fulongfei/analysis/hiv/generateCons_test/check_same_as_ref.pl $file $tvc_freq_file $outfile`;
}