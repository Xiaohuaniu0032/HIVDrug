use strict;
use warnings;
use File::Basename;
use FindBin qw/$Bin/;

# /data/langna/projects/0.RD/hivtools_v2/hivtest_v1/output/v0/IonXpress_002/consensus

my $ref = "/data/fulongfei/git_repo/HIVDrug/ref/K03455.fasta";
my $resdir = "/data/langna/projects/0.RD/hivtools_v2/hivtest_v1/output/v0";

my ($outdir) = @ARGV;

my @summary = glob "$resdir/*/consensus/summary.json";

for my $val (@summary){
	#print "$val\n";
	my $dir = dirname($val);
	my $freq_1 = "$dir/Freq.0.01.consensus.fasta";
	my $freq_2 = "$dir/Freq.0.02.consensus.fasta";
	my $freq_5 = "$dir/Freq.0.05.consensus.fasta";
	my $freq_10 = "$dir/Freq.0.1.consensus.fasta";
	my $freq_20 = "$dir/Freq.0.2.consensus.fasta";

	my @files;
	push @files, $freq_1;
	push @files, $freq_2;
	push @files, $freq_5;
	push @files, $freq_10;
	push @files, $freq_20;

	my $dd = dirname($dir); # /data/langna/projects/0.RD/hivtools_v2/hivtest_v1/output/v0/IonXpress_002
	my $name = basename($dd); # IonXpress_002
	my $full_name = $name."_rawlib";
	print "$full_name\n";
	
	if (!-d "$outdir/$full_name/nextalign"){
		`mkdir -p $outdir/$full_name/nextalign`;
	}

	my $runsh = "$outdir/$full_name/nextalign/run.sh";
	my $input_fa = "$outdir/$full_name/nextalign/input.fa";

	open CONS, ">$input_fa" or die;
	for my $file (@files){
		open IN, "$file" or die;
		my $h = <IN>;
		my $seq = <IN>;
		print CONS "$h";
		print CONS "$seq\n";
	}
	close CONS;

	open O, ">$runsh" or die;
	#print O "cat $freq_1 $freq_2 $freq_5 $freq_10 $freq_20 >$input_fa\n";
	print O "/data/fulongfei/git_repo/ClustalO_2019nCoV/bin/nextalign --sequences\=$input_fa --reference\=$ref --output-dir\=$outdir/$full_name/nextalign --output-basename\=$full_name --include-reference --in-order\n";
	close O;
}
