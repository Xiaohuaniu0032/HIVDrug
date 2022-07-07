use strict;
use warnings;
use Data::Dumper;
use FindBin qw/$Bin/;

my ($align_file,$outfile) = @ARGV;

open O, ">$outfile" or die;

my @seq_name;
my %seq_info;
open IN, "$align_file" or die;
while (<IN>){
	chomp;
	# header
	my $header = (split /\s+/, $_)[0];
	$header =~ s/^\>//;
	my $seq = <IN>;
	chomp $seq;
	push @seq_name, $header;
	$seq_info{$header} = $seq;
}
close IN;

#print "@seq_name\n";

print O "Chr\tPos\t";

for my $seq (@seq_name){
	print O "\t$seq";
}
print O "\n";


my $ref_name = $seq_name[0];
my $ref_fa = $seq_info{$ref_name};
my $first_query_name = $seq_name[1];
my $first_query_seq  = $seq_info{$first_query_name};

my $len_ref = length($ref_fa);
my $len_query = length($first_query_seq);

my @ref_fa = split //, $ref_fa;
my @query_fa = split //, $first_query_seq;

# determine first and last ATCG base in query seq
my @pos;
my $pos_idx = 0;
for my $base (@query_fa){
	$pos_idx += 1;
	if ($base =~ /[ATCG]/){
		push @pos, $pos_idx;
	}
}

#print "@pos\n";

my $start_pos = $pos[0];
my $end_pos = $pos[-1];

print "$start_pos\t$end_pos\n";

my %seq_NNN;
for my $name (@seq_name){
	my $seq = $seq_info{$name};
	my @seq = split //, $seq;
	$seq_NNN{$name} = \@seq;
}

for (my $i=0;$i<=$len_ref-1;$i++){
	my @nt;
	for my $name (@seq_name){
		my @base = @{$seq_NNN{$name}};
		my $base = $base[$i];
		push @nt, $base;
	}
	#print "@nt\n";

	my $ref_base = shift @nt;

	my $first_base = $nt[0];
	my $same_flag = 1;
	for my $base (@nt){
		if ($base ne $first_base){
			$same_flag = 0;
		}
	}

	my $if_same;
	if ($same_flag == 1){
		# all same
		$if_same = "Same";
	}else{
		$if_same = "NotSame";
	}

	my $idx = $i + 1;
	next if ($idx < $start_pos || $idx > $end_pos);
	my $base_str = join("\t",@nt);
	print O "K03455.1\t$idx\t$ref_base\t$base_str\t$if_same\n";
	#print "$idx\n";
}
close O;
