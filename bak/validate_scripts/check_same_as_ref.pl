use strict;
use warnings;

my ($parsed_file,$tvc_igv_freq_file,$outfile) = @ARGV;

open O, ">$outfile" or die;


my %tvc_var_info;
open VAR, "$tvc_igv_freq_file" or die;
while (<VAR>){
	chomp;
	my @arr = split /\t/;
	next if ($arr[-1] eq "NA");
	my $pos = $arr[1];
	my $ref = $arr[2];
	my $alt = $arr[3];
	my $igv_freq = $arr[-2];

	my $tvc_freq_tmp = $arr[-1];
	my $tvc_freq;
	if ($tvc_freq_tmp =~ /\;/){
		my @tmp = split /\;/, $tvc_freq_tmp;
		for my $v (@tmp){
			$tvc_freq += $v;
		}
	}else{
		$tvc_freq = $tvc_freq_tmp;
	}

	my $var = "$pos\:$ref\:$alt\:$tvc_freq";
	push @{$tvc_var_info{$pos}}, $var; # one pos may have 1/2/... vars
}
close VAR;

open IN, "$parsed_file" or die;
my $header = <IN>;
chomp $header;
my @h = split /\t/, $header;
pop @h;
my $h = join("\t",@h);
print O "$h\t";
print O "SameAsRef\tTVC_var_info\n";


while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $pos = $arr[1];
	
	my $ref_base = $arr[2];
	
	my $freq2_base = $arr[3];
	my $freq5_base = $arr[4];
	my $freq10_base = $arr[5];
	my $freq15_base = $arr[6];
	my $freq20_base = $arr[7];

	my $same_as_ref;
	if ($freq2_base ne $ref_base || $freq5_base ne $ref_base || $freq10_base ne $ref_base || $freq15_base ne $ref_base || $freq20_base ne $ref_base){
		# 存在突变
		$same_as_ref = "NotSameAsRef";
	}else{
		$same_as_ref = "SameAsRef";
	}

	my $tvc_var;
	if (exists $tvc_var_info{$pos}){
		# this pos has a tvc call
		my $tvc_var_aref = $tvc_var_info{$pos};
		$tvc_var = join(";", @{$tvc_var_aref});
	}else{
		$tvc_var = "NA";
	}

	print O "$arr[0]\t$pos\t$ref_base\t$freq2_base\t$freq5_base\t$freq10_base\t$freq15_base\t$freq20_base\t$same_as_ref\t$tvc_var\n";
}
close IN;
close O;