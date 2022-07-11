use strict;
use warnings;
use Data::Dumper;
use FindBin qw/$Bin/;

my ($bam,$name,$ref_fa,$bed,$outdir) = @ARGV;

my $runsh = "$outdir/run\_$name\.sh";
open O, ">$runsh" or die;

print O "# pileup2vaf\n";
print O "/usr/bin/python3 /data/fulongfei/git_repo/pileup2vaf/pileup2vaf.py -bam $bam -fa $ref_fa -bed $bed -od $outdir -n $name\n";

print O "# "

