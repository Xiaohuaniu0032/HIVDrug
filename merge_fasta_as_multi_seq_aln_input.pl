use strict;
use warnings;
use Data::Dumper;
use FindBin qw/$Bin/;
use File::Basename;

#my $resdir = "/data/fulongfei/analysis/hiv/JSCDC";
my $ref_fa = "/data/fulongfei/git_repo/HIVDrug/ref/K03455.fasta";

my ($indir,$name) = @ARGV;

my @runsh = glob "$indir/$name/*.sh";

for my $sh (@runsh){
	my $dir = dirname($sh);
	my $name = basename($dir);

	my $input_fa = "$dir/nextalign/$name\.aln.input.fa";

	open O, ">$input_fa" or die;
	
	# write ref fasta
	print O ">K03455.1\n";
	my @ref_fa;
	open REF, "$ref_fa" or die;
	while (<REF>){
		chomp;
		next if /^$/;
		next if /^\>/;
		push @ref_fa, $_;
	}
	close REF;
	my $ref_fa = join("",@ref_fa);
	print O "$ref_fa\n";

	# write CE fasta
	my @CE_file = glob "$dir/CE/*.fas";
	my $CE_file = $CE_file[0];

	`/usr/bin/dos2unix $CE_file`;

	my @CE_fa;
	print O "\>$name\_CE\n";
	open CE, "$CE_file" or die;
	while (<CE>){
		chomp;
		next if /^$/;
		next if /^\>/;
		push @CE_fa, $_;
	}
	close CE;
	my $ce_fa = join("",@CE_fa);
	print O "$ce_fa\n";

	# wirte IRMA
	my $irma_file = "$dir/IRMA/Boxin_Region/consensusFreq20/amended_consensus/consensusFreq20.fa";
	open IRMA, "$irma_file" or die;
	while (<IRMA>){
		chomp;
		if (/^\>/){
			print O "\>$name\_consensusFreq20\_IRMA\n";
		}else{
			print O "$_\n";
		}
	}
	close IRMA;
	
	# write freebayesConsensus_v2
	#my @freq = qw/2 5 10 15 20/;
	my @freq = qw/20/;
	for my $freq (@freq){
		my $Freq_STR = "Freq".$freq; # Freq20
		my $freebayes_fa = "$dir/freebayesConsensus_v2/$Freq_STR/$name\.freebayesConsensus\.$Freq_STR\.fasta";

		open FA, "$freebayes_fa" or die;
		while (<FA>){
			chomp;
			print O "$_\n";
		}
		close FA;
	}

	close O;
}
