use strict;
use warnings;

my ($allele_file,$pileup_file) = @ARGV;

my %allele_var;
open IN, "$allele_file" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $pos = $arr[1];
	my $ref = $arr[2];
	my $alt = $arr[3];
	my $freq = $arr[6];

	# alleles.xls文件中，有些位点会有重复
	# HIV-1   2806    T       C       Heterozygous    -       38.5    5864.9  -       SNP     Novel
	# HIV-1   2806    TAA     CAG     Heterozygous    -       60.9    5864.9  -       MNP     Novel

	my $type = $arr[9]; # MNP/SNP
	if ($type eq "SNP"){
		my $var = "$arr[1]\t$arr[2]\t$arr[3]"; # pos/ref/alt
		push @{$allele_var{$var}}, $freq;
	}elsif ($type eq "MNP"){
		my $len = length($ref);
		my @ref_arr = split //, $ref;
		my @alt_arr = split //, $alt;
		for my $i (0..$len-1){
			my $ref_base = $ref_arr[$i];
			my $alt_base = $alt_arr[$i];
			if ($ref_base ne $alt_base){
				my $rel_pos = $pos + $i;
				my $var = "$rel_pos\t$ref_base\t$alt_base";
				push @{$allele_var{$var}}, $freq;
			}
		}
	}else{
		next; # skip ins/del/complex type
	}
}
close IN;
#print(Dumper(\%allele_var));


my %pileup_var;
open IN, "$pileup_file" or die;
my $h = <IN>;
chomp $h;
print "$h\tvariantCaller_Freq\n";
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $ref = $arr[2];
	my $alt = $arr[3];
	my $ref_len = length($ref);
	my $alt_len = length($alt);
	next if ($ref_len != $alt_len); # skip indel
	if ($ref_len == 1 and $alt_len == 1){
		my $var = "$arr[1]\t$ref\t$alt"; # pos/ref/alt
		if (exists $allele_var{$var}){
			my @tvc_freq = @{$allele_var{$var}};
			my $tvc_freq = join(";",@tvc_freq);
			print "$_\t$tvc_freq\n";
		}else{
			print "$_\tNA\n"; # if do not have tvc call, then this pos's tvc freq will be NA
		}
	}
}
close IN;

