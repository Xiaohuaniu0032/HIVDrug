use strict;
use warnings;
use File::Basename;

# IRAM	/data/langna/projects/0.RD/hivtools_v2/hivtools_v1/output/IonXpress_002/IRMA/consensusFreq20/amended_consensus
# generateCons	/data/langna/projects/0.RD/hivtools_v2/hivtest_v1/output/v0/IonXpress_002/consensus/Freq.0.2.consensus.fasta

my $ref = "/data/fulongfei/git_repo/HIVDrug/ref/K03455.fasta";

my ($outdir) = @ARGV;

my $irma_dir = "/data/langna/projects/0.RD/hivtools_v2/hivtools_v1/output";
my @IRMA_fa = glob "$irma_dir/*/IRMA/consensusFreq20/amended_consensus/consensusFreq20.fa";

my $cons_dir = "/data/langna/projects/0.RD/hivtools_v2/hivtest_v1/output/v0";

for my $fa (@IRMA_fa){
	#print "$fa\n";
	my $d = dirname($fa); # /data/langna/projects/0.RD/hivtools_v2/hivtools_v1/output/IonXpress_002/IRMA/consensusFreq20/amended_consensus
	my $dd = dirname($d); # /data/langna/projects/0.RD/hivtools_v2/hivtools_v1/output/IonXpress_002/IRMA/consensusFreq20
	my $ddd = dirname($dd); # /data/langna/projects/0.RD/hivtools_v2/hivtools_v1/output/IonXpress_002/IRMA
	my $dddd = dirname($ddd); # /data/langna/projects/0.RD/hivtools_v2/hivtools_v1/output/IonXpress_002
	my $name = basename($dddd); # IonXpress_002

	if (!-d "$outdir/$name"){
		`mkdir $outdir/$name`;
	}

	my $cons_fa = "$cons_dir/$name/consensus/Freq.0.2.consensus.fasta";
	

	my $input_fa = "$outdir/$name/input.fa";
	
	my $runsh = "$outdir/$name/run.sh";
	open O, ">$runsh" or die;
	print O "cat $fa $cons_fa >$input_fa\n";
	print O "/data/fulongfei/git_repo/ClustalO_2019nCoV/bin/nextalign --sequences\=$input_fa --reference\=$ref --output-dir\=$outdir/$name --output-basename\=$name --include-reference --in-order\n";
	my $aln_fa = "$outdir/$name/$name\.aligned.fasta";
	my $parsed_file = "$outdir/$name/$name\.parsed.txt";
	print O "perl /data/fulongfei/git_repo/HIVDrug/scripts/parse_nextalign.pl $aln_fa $parsed_file\n";
	close O;
}
