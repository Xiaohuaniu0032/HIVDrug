use strict;
use warnings;
use File::Basename;


my ($pileup_dir) = @ARGV;

my @var = glob "$pileup_dir/*/pileup2vaf/*.variants.xls";

my $tvc_dir = $pileup_dir;

for my $file (@var){
	my $d = dirname($file); # /data/fulongfei/analysis/hiv/variantCaller/IonXpress_001_rawlib/pileup2vaf
	my $dd = dirname($d); # /data/fulongfei/analysis/hiv/variantCaller/IonXpress_001_rawlib
	my $name = basename($dd); # IonXpress_001_rawlib
	print "$name\n";

	my $allele = "$tvc_dir/$name/variantCaller/alleles.xls";
	
	my $runsh = "$tvc_dir/$name/freq_cmp.sh";
	open O, ">$runsh" or die;
	print O "perl /data/fulongfei/analysis/hiv/final_tvc_json_test/compare_variantCaller_and_pileup2vaf.pl $allele $file >$tvc_dir/$name/freq_cmp.txt\n";
	close O;
}

