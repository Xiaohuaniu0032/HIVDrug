use strict;
use warnings;
use File::Basename;
use FindBin qw/$Bin/;

my ($outdir) = @ARGV;

my $ref = "$Bin/ref/K03455.fasta";
my $samtools = "/usr/bin/samtools";

my %codon_aa;
my $codon_file = "$Bin/codon.list";

open IN, "$codon_file" or die;
while (<IN>){
	chomp;
	next if /^\#/;
	next if /^$/;
	next if /^Amino/; # skip header
	# Isoleucine      I       ATT, ATC, ATA
	my @arr = split /\t/;
	my $aa  = "$arr[1]\t$arr[0]";
	my $dna = $arr[2];
	my @dna = split /\,/, $dna;
	for my $dna (@dna){
		$dna =~ s/^\s+//;
		$codon_aa{$dna} = $aa; # AAT => "I \t Isoleucine"
	}
}
close IN;


my %AA_Pos;
$AA_Pos{"Protease"} = "2253,2549";
$AA_Pos{"RT"}       = "2550,4229";
$AA_Pos{"Integrase"} = "4230,5096";

my @AA = qw/Protease RT Integrase/;

for my $aa (@AA){
	my $pos = $AA_Pos{$aa};
	my @pos = split /\,/, $pos;
	my $sp = $pos[0];
	my $ep = $pos[1];

	my $outfile = "$outdir/$aa\.$sp\.$ep\.fasta";
	my $cmd = "$samtools faidx $ref K03455.1\:$sp\-$ep >$outfile";
	system($cmd);

	my $seq = cat_fa($outfile);
	#print "$seq\n";
	my $len = length($seq);

	my $idx = 0;
	for (my $i=0; $i<=$len-3; $i+=3){
		$idx += 1;
		my $codon = substr($seq,$i,3);
		
		my $val;
		if (exists $codon_aa{$codon}){
			$val = $codon_aa{$codon};
		}else{
			$val = "NA";
		}
		
		print "$aa\t$idx\t$codon\t$val\n";
	}

}


sub cat_fa{
	my ($fa) = @_;
	my @seq;
	open IN, "$fa" or die;
	while (<IN>){
		chomp;
		next if /^\>/;
		next if /^$/;
		push @seq, $_;
	}
	close IN;
	my $seq = join("",@seq);

	return($seq);
}

# HIVdb FAQ
# https://hivdb.stanford.edu/pages/FAQ/FAQ.html

# mut-scores
# https://hivdb.stanford.edu/dr-summary/mut-scores/PI/
# https://hivdb.stanford.edu/dr-summary/mut-scores/NRTI/
# https://hivdb.stanford.edu/dr-summary/mut-scores/NRTI/
# https://hivdb.stanford.edu/dr-summary/mut-scores/INSTI/


# Ref is K03455.1
# CDS [2358..5096] [pol polyprotein] [codon_start=1]

# Protease                  [2253,2549]
# Reverse Transcriptase     [2550,4229]
# Integrase                 [4233,5096]


# https://hivdb.stanford.edu/hivdb/by-mutations/
# https://www.hiv.lanl.gov/content/sequence/HIV/REVIEWS/HXB2.html
# https://www.ncbi.nlm.nih.gov/nuccore/K03455.1

# Protease
# >>>Pol p10 Protease
# >>>PQVTLWQRPL VTIKIGGQLK EALLDTGADD TVLEEMSLPG RWKPKMIGGI GGFIKVRQYD QILIEICGHK AIGTVLVGPT PVNIIGRNLL TQIGCTLNF    99

# Reverse Transcriptase
# >>>Pol p66 Reverse Transcriptase (RT/RNAse)
# >>>PISPIETVPV KLKPGMDGPK VKQWPLTEEK IKALVEICTE MEKEGKISKI GPENPYNTPV FAIKKKDSTK WRKLVDFREL NKRTQDFWEV QLGIPHPAGL  100
# >>>KKKKSVTVLD VGDAYFSVPL DEDFRKYTAF TIPSINNETP GIRYQYNVLP QGWKGSPAIF QSSMTKILEP FRKQNPDIVI YQYMDDLYVG SDLEIGQHRT  200
# >>>KIEELRQHLL RWGLTTPDKK HQKEPPFLWM GYELHPDKWT VQPIVLPEKD SWTVNDIQKL VGKLNWASQI YPGIKVRQLC KLLRGTKALT EVIPLTEEAE  300
# >>>LELAENREIL KEPVHGVYYD PSKDLIAEIQ KQGQGQWTYQ IYQEPFKNLK TGKYARMRGA HTNDVKQLTE AVQKITTESI VIWGKTPKFK LPIQKETWET  400
# >>>WWTEYWQATW IPEWEFVNTP PLVKLWYQLE KEPIVGAETF YVDGAANRET KLGKAGYVTN RGRQKVVTLT DTTNQKTELQ AIYLALQDSG LEVNIVTDSQ  500
# >>>YALGIIQAQP DQSESELVNQ IIEQLIKKEK VYLAWVPAHK GIGGNEQVDK LVSAGIRKVL                                              560

# Integrase
# >>>Pol p31 Integrase
# >>>FLDGIDKAQD EHEKYHSNWR AMASDFNLPP VVAKEIVASC DKCQLKGEAM HGQVDCSPGI WQLDCTHLEG KVILVAVHVA SGYIEAEVIP AETGQETAYF  100
# >>>LLKLAGRWPV KTIHTDNGSN FTGATVRAAC WWAGIKQEFG IPYNPQSQGV VESMNKELKK IIGQVRDQAE HLKTAVQMAV FIHNFKRKGG IGGYSAGERI  200
# >>>VDIIATDIQT KELQKQITKI QNFRVYYRDS RNPLWKGPAK LLWKGEGAVV IQDNSDIKVV PRRKAKIIRD YGKQMAGDDC VASRQDED               288


