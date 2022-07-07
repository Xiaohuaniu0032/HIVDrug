use strict;
use warnings;
use File::Basename;

my ($bam_dir,$outdir) = @ARGV;

my $bed = '/data/fulongfei/analysis/hiv/re_analysis_new_ref/POL.bed';
my $ref = '/data/fulongfei/git_repo/HIVDrug/ref/K03455.fasta';

my @bam = glob "$bam_dir/*.bam";
for my $bam (@bam){
	my $name = (split /\./, basename($bam))[0];
	print "$name\n";
	
	if (!-d "$outdir/$name"){
		`mkdir $outdir/$name`;
	}

	if (!-d "$outdir/$name/pileup2vaf"){
		`mkdir $outdir/$name/pileup2vaf`;
	}

	my $pileup_sh = "$outdir/$name/pileup2vaf/run.sh";
	open O, ">$pileup_sh" or die;
	print O "python3 /data/fulongfei/git_repo/pileup2vaf/pileup2vaf.py -bam $bam -fa $ref -bed $bed -od $outdir/$name/pileup2vaf -n $name\n";
	#print O "perl /data/fulongfei/git_repo/pileup2vaf/scripts/filter_var_by_freq.pl $outdir/$name/pileup2vaf/$name\.variants.xls 0.04 $outdir/$name/pileup2vaf/$name\.variants.filter.xls\n";
	close O;
}
