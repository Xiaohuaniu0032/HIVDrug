use strict;
use warnings;
use File::Basename;
use FindBin qw/$Bin/;

# inpput file
# 1) drug.aa.list.xls
# Protease        1       CCT     P       Proline
# Protease        2       CAG     Q       Glutamine
# Protease        3       GTC     V       Valine


# 2) db/ScoresPI_1655774816824.tsv
# Rule    Position        AA      BIC     CAB     DTG     EVG     RAL
# H51Y    51      Y       10      15      10      15      15
# T66A    66      A       0       0       0       60      15

# Temp Format:Integrase H51Y,H51K 1-3 CAT TAT,TAC

# output format
# Chr/Pos/ID/Ref/Alt
# K03455.1/1/Integrase,H51Y,H51K/C/T
# K03455.1/3/Integrase,H51Y,H51K/T/C

my ($score_tsv_file,)
