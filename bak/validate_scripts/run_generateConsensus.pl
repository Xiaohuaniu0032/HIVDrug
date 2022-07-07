use strict;
use warnings;
use File::Basename;

my ($tvc_dir,$outdir) = @ARGV;

# /data/langna/projects/0.RD/hivtools_v2/hivtest_v1/output/v0/IonXpress_002/TVC/TSVC_variants.genome.vcf
#my @gvcf = glob "$tvc_dir/*/variantCaller/TSVC_variants.genome.vcf";
my @gvcf = glob "$tvc_dir/*/TVC/TSVC_variants.genome.vcf";

for my $gvcf (@gvcf){
	my $d = dirname($gvcf); # /data/langna/projects/0.RD/hivtools_v2/hivtest_v1/output/v0/IonXpress_002/TVC
	my $dd = dirname($d); # /data/langna/projects/0.RD/hivtools_v2/hivtest_v1/output/v0/IonXpress_002
	my $name = basename($dd); # IonXpress_002
	
	my @freq = qw/0.02 0.05 0.1 0.15 0.2/;
	for my $freq (@freq){
		my $freq_tmp = $freq * 100;
		my $dir_tmp = "$outdir/$name/generateConsensus/Freq".$freq_tmp;
		if (!-d $dir_tmp){
			`mkdir -p $dir_tmp`;
		}

		my $runsh = "$dir_tmp/run.sh";
		open O, ">$runsh" or die;
		my $fa_h = $name."_generateCons_Freq".$freq_tmp;
		my $cons_fa = "$dir_tmp/$name\.consensus\.Freq"."$freq_tmp\.fasta";
		print O "python /data/fulongfei/git_repo/HIVDrug/generateConsensus/gvcf_to_fasta_v2.py --gvcf-file $gvcf --input-fasta-file /data/fulongfei/git_repo/HIVDrug/ref/K03455.fasta --region-bed-file /data/fulongfei/git_repo/HIVDrug/BED/new_BED/WG00471_HIV_pol.20201027_Designed.NewRef.bed --output-fasta-file $cons_fa --alias-contig $fa_h --min-dp 20 --process-contig K03455.1 --major-allele-only 0 --min-hpindel-var-freq $freq --min-non-hpindel-var-freq $freq\n";
		close O;

		`chmod 755 $runsh`;
		
		print "Runing $runsh\.\.\.\n";
		system($runsh);
	}
}

