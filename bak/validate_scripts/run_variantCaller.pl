use strict;
use warnings;
use File::Basename;

my ($bam_dir,$outdir) = @ARGV;

my $bed = '/data/fulongfei/analysis/hiv/re_analysis_new_ref/POL.bed';
my $hotspot_vcf = "/data/fulongfei/git_repo/HIVDrug/scripts/make_hs_vcf_from_bed/POL.hotspot.vcf";
my $ref = '/data/fulongfei/git_repo/HIVDrug/ref/K03455.fasta';
#my $config = '/data/fulongfei/tools/variantCaller/pluginMedia/parameter_sets/ampliseq_somatic_lowstringency_530_parameters.json';
#my $config = '/data/fulongfei/analysis/hiv/re_analysis_new_ref/ampliseq_somatic_lowstringency_530_parameters.json'; # snp limit = 100 and PRESHIT=0.3
my $config = "/data/fulongfei/analysis/hiv/final_tvc_json_test/ampliseq_somatic_lowstringency_530_parameters.json";
my @bam = glob "$bam_dir/*.bam";
for my $bam (@bam){
	my $name = (split /\./, basename($bam))[0];
	print "$name\n";
	
	if (!-d "$outdir/$name"){
		`mkdir $outdir/$name`;
	}

	if (!-d "$outdir/$name/variantCaller"){
		`mkdir $outdir/$name/variantCaller`;
	}

	my $tvc_sh = "$outdir/$name/variantCaller/run.sh";
	open O, ">$tvc_sh" or die;
	print O "python /data/fulongfei/tools/variantCaller/bin/variant_caller_pipeline.py -b $bed -i $bam -r $ref -o $outdir/$name/variantCaller -p $config -m /data/fulongfei/tools/variantCaller/share/TVC/sse/motifset.txt --generate-gvcf on --hotspot-vcf $hotspot_vcf\n";
	my $TSVC_variants = "$outdir/$name/variantCaller/TSVC_variants.vcf";
	my $alleles_xls   = "$outdir/$name/variantCaller/alleles.xls";
	my $variants_xls  = "$outdir/$name/variantCaller/variants.xls";

	print O "python /data/fulongfei/tools/variantCaller/scripts/generate_variant_tables.py --suppress-no-calls on --input-vcf $TSVC_variants --region-bed $bed --alleles2-xls $alleles_xls --output-xls $variants_xls\n";
	close O;
	`chmod 755 $tvc_sh`;
}
