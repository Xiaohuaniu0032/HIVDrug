use strict;
use warnings;

my ($bed) = @ARGV;
my $shift_pos = 454;

# WG00471_HIV_pol.20201027_Designed.bed
my $ref_name = "K03455.1";

open IN, "$bed" or die;
my $h = <IN>;
print "$h";
while (<IN>){
	chomp;
	my @arr = split /\t/, $_;
	my $sp = $arr[1];
	my $ep = $arr[2];
	my $len = $ep - $sp;

	my $new_sp = $sp + $shift_pos;
	my $new_ep = $new_sp + $len;
	print "$ref_name\t$new_sp\t$new_ep\t$arr[3]\t$arr[4]\t$arr[5]\n";
}
close IN;
