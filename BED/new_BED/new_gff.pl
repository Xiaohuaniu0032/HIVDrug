use strict;
use warnings;

my ($gff) = @ARGV;
# genes.gff

my $shift_pos = 454;
my $ref_name = "K03455.1";


open IN, "$gff" or die;
while (<IN>){
	chomp;
	if (/^\#/){
		print "$_\n";
	}else{
		my @arr = split /\t/, $_;
		my $sp = $arr[3];
		my $ep = $arr[4];
		my $len = $ep - $sp;

		my $new_sp = $sp + $shift_pos;
		my $new_ep = $new_sp + $len;
		print "$ref_name\t$arr[1]\t$arr[2]\t$new_sp\t$new_ep\t$arr[5]\t$arr[6]\t$arr[7]\t$arr[8]\n";
	}
}
close IN;
