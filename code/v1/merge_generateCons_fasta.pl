use strict;
use warnings;
use File::Basename;

my ($resdir,$name,$nextalign_dir) = @ARGV;

my $cons_freq2 = "$resdir/Freq2/$name\.consensus.Freq2.fasta";
my $cons_freq5 = "$resdir/Freq5/$name\.consensus.Freq5.fasta";
my $cons_freq10 = "$resdir/Freq10/$name\.consensus.Freq10.fasta";
my $cons_freq15 = "$resdir/Freq15/$name\.consensus.Freq15.fasta";
my $cons_freq20 = "$resdir/Freq20/$name\.consensus.Freq20.fasta";

my @fa;
push @fa, $cons_freq2;
push @fa, $cons_freq5;
push @fa, $cons_freq10;
push @fa, $cons_freq15;
push @fa, $cons_freq20;

my $input_fa = "$nextalign_dir/input.fa";
open FA, ">$input_fa" or die "$!\n";
for my $fa_file (@fa){
	open FASTA, "$fa_file" or die "$!\n";
	my $h = <FASTA>;
	chomp $h;
	my $seq = <FASTA>;
	chomp $seq;
	print FA "$h\n$seq\n";
	close FASTA;
}
close FA;