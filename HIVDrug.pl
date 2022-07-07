use strict;
use warnings;
use Data::Dumper;
use FindBin qw/$Bin/;

my ($gvcf,$name,$ref_fa,$bed,$outdir) = @ARGV;

my $runsh = "$outdir/run.sh";
open O, ">$runsh" or die;
print O "python $Bin/scripts/gvcf_to_fasta.py --gvcf-file $gvcf --input-fasta-file $ref_fa --region-bed-file $bed --output-fasta-file $outdir/$name\.consensus.fasta --alias-contig $name --min-dp 20 --process-contig K03455.1 --major-allele-only 0 --min-hpindel-var-freq 0.6 --min-non-hpindel-var-freq 0.02\n";
my $cons_fa = "$outdir/$name\.consensus.fasta";
print O "/data/fulongfei/git_repo/ClustalO_2019nCoV/bin/nextalign --sequences\=$cons_fa --reference\=$ref_fa --output-dir\=$outdir --output-basename\=$name --include-reference\n";
close O;
