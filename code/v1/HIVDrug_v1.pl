use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw/$Bin/;
use List::Util qw(sum);
use POSIX qw(ceil);
use POSIX qw(floor);

# variantCaller 5.12
# generateConsensus
# IRMA
# pileup2vaf
# nextalign
# freebayes

my ($bam,$name,$outdir) = @ARGV;

my $python2 = "/usr/bin/python";
my $python3 = "/usr/bin/python3";

my $unmerged_bed = "/data/fulongfei/git_repo/HIVDrug/BED/new_BED/WG00471_HIV_pol.20201027_Designed.NewRef.bed";
my $merged_bed = "/data/fulongfei/git_repo/HIVDrug/BED/new_BED/POL.bed";
my $ref = "/data/fulongfei/git_repo/HIVDrug/ref/K03455.fasta";
my $tvcJson = "/data/fulongfei/analysis/hiv/final_tvc_json_test/ampliseq_somatic_lowstringency_530_parameters.json";
my $tvcSse = "/data/fulongfei/tools/variantCaller/share/TVC/sse/motifset.txt";
my $hotspot_vcf = "/data/fulongfei/git_repo/HIVDrug/scripts/make_hs_vcf_from_bed/POL.hotspot.vcf";

my $freebayes = "/data/fulongfei/tools/freebayes-1.3.4-linux-static-AMD64";

if (!-d "$outdir/$name/freebayes"){
	`mkdir -p $outdir/$name/freebayes`;
}

if (!-d "$outdir/$name/variantCaller"){
	`mkdir -p $outdir/$name/variantCaller`;
}

if (!-d "$outdir/$name/pileup2vaf"){
	`mkdir -p $outdir/$name/pileup2vaf`;
}

if (!-d "$outdir/$name/generateConsensus"){
	`mkdir -p $outdir/$name/generateConsensus`;
}

if (!-d "$outdir/$name/nextalign"){
	`mkdir -p $outdir/$name/nextalign`;
}

if (!-d "$outdir/$name/IRMA"){
	`mkdir -p $outdir/$name/IRMA`;
}

if (!-d "$outdir/$name/freebayesConsensus"){
	`mkdir -p $outdir/$name/freebayesConsensus`;
}

if (!-d "$outdir/$name/freebayesConsensus_v2"){
	`mkdir -p $outdir/$name/freebayesConsensus_v2`;
}


my $runsh = "$outdir/$name/run\_$name.sh";

open O, ">$runsh" or die;
print O "\# variantCaller 5.12\n";
print O "\# $python2 /data/fulongfei/tools/variantCaller/bin/variant_caller_pipeline.py --region-bed $merged_bed --primer-trim-bed $unmerged_bed --input-bam $bam --reference-fasta $ref --output-dir $outdir/$name/variantCaller --parameters-file $tvcJson -m $tvcSse --generate-gvcf on --hotspot-vcf $hotspot_vcf\n\n";
my $TSVC_variants_vcf = "$outdir/$name/variantCaller/TSVC_variants.vcf";
my $alleles_xls = "$outdir/$name/variantCaller/alleles.xls";
my $variants_xls = "$outdir/$name/variantCaller/variants.xls";
print O "\# $python2 /data/fulongfei/tools/variantCaller/scripts/generate_variant_tables.py --suppress-no-calls on --input-vcf $TSVC_variants_vcf --region-bed $unmerged_bed --alleles2-xls $alleles_xls --output-xls $variants_xls\n\n";


print O "\# freebayes\n";
my @ploidy = qw/2 3/;
for my $p (@ploidy){
	next if ($p == 3);
	my $pp = "ploidy".$p;
	#print O "$freebayes --bam $bam --fasta-reference $ref -t $merged_bed --vcf $outdir/$name/freebayes/$name\.freebayes.$pp\.vcf -F 0.02 -C 1 --ploidy $p --pooled-continuous --gvcf $outdir/$name/freebayes/$name\.freebayes.$pp\.gvcf\n\n";
	print O "\# $freebayes --bam $bam --fasta-reference $ref -t $merged_bed -F 0.02 -C 1 --limit-coverage 3000 --ploidy $p --pooled-continuous --gvcf \>$outdir/$name/freebayes/$name\.freebayes\.$pp\.gvcf\n\n";
}

print O "\# freebayesConsensus\n";
my @freq = qw/0.02 0.05 0.1 0.15 0.2/;
my $freebayes_gvcf = "$outdir/$name/freebayes/$name\.freebayes.ploidy2.gvcf";
for my $f (@freq){
	my $freq_int = $f * 100;
	my $freq_new = "Freq".$freq_int;
	if (!-d "$outdir/$name/freebayesConsensus/$freq_new"){
		`mkdir -p $outdir/$name/freebayesConsensus/$freq_new`;
	}
	print O "\# perl /data/fulongfei/git_repo/HIVDrug/freebayes_consensus_fa.pl -gvcf $freebayes_gvcf -n $name -ref $ref -bed $merged_bed -d 20 -snp_f $f -indel_f 0.6 -outdir $outdir/$name/freebayesConsensus/$freq_new\n\n";
}

print O "\# freebayesConsensus_v2\n";
for my $f (@freq){
	my $freq_int = $f * 100;
	my $freq_new = "Freq".$freq_int;
	if (!-d "$outdir/$name/freebayesConsensus_v2/$freq_new"){
		`mkdir -p $outdir/$name/freebayesConsensus_v2/$freq_new`;
	}
	print O "perl /data/fulongfei/git_repo/HIVDrug/freebayes_consensus_fa_v2.pl -gvcf $freebayes_gvcf -n $name -ref $ref -bed $merged_bed -d 20 -snp_f $f -indel_f 0.6 -outdir $outdir/$name/freebayesConsensus_v2/$freq_new\n\n";
}


print O "\# pileup2vaf\n";
print O "\# $python3 /data/fulongfei/git_repo/pileup2vaf/pileup2vaf.py -bam $bam -fa $ref -bed /data/fulongfei/analysis/hiv/re_analysis_new_ref/POL.bed -od $outdir/$name/pileup2vaf -n $name\n\n";

my $gvcf = "$outdir/$name/variantCaller/TSVC_variants.genome.vcf";

print O "\# generateConsensus\n";
undef @freq;
@freq = qw/0.2/;
for my $freq (@freq){
	my $freq_tmp = $freq * 100;
	
	my $dir_tmp = "$outdir/$name/generateConsensus/Freq".$freq_tmp; # Freq2/Freq5/Freq10/Freq15/Freq20
	if (!-d $dir_tmp){
		`mkdir -p $dir_tmp`;
	}

	my $fa_h = $name."_generateCons_Freq".$freq_tmp;
	my $cons_fa = "$dir_tmp/$name\.consensus\.Freq"."$freq_tmp\.fasta";
	print O "\# $python2 /data/fulongfei/git_repo/HIVDrug/generateConsensus/gvcf_to_fasta_v2.py --gvcf-file $gvcf --input-fasta-file $ref --region-bed-file $unmerged_bed --output-fasta-file $cons_fa --alias-contig $fa_h --min-dp 20 --process-contig K03455.1 --major-allele-only 0 --min-hpindel-var-freq $freq --min-non-hpindel-var-freq $freq\n\n";
}


print O "\# IRMA\n";
my $fq = "$outdir/$name/IRMA/$name\.fastq";
print O "\# samtools fastq $bam >$fq\n\n";

# Boxin region
if (!-d "$outdir/$name/IRMA/Boxin_Region"){
	`mkdir -p $outdir/$name/IRMA/Boxin_Region`;
}
for my $freq (@freq){
	my $freq_tmp = $freq * 100;
	my $freq_str = "consensusFreq".$freq_tmp;
	my $cons_dir = "$outdir/$name/IRMA/Boxin_Region/$freq_str";
	my $IRMA_module_dir = "HIV1freq".$freq_tmp."_Boxin-pgm-hq"; # HIV1freq20_Boxin
	print O "\# /data/fulongfei/git_repo/flu-amd/IRMA $IRMA_module_dir $fq $cons_dir\n\n";
}
# JSCDC_Segment
if (!-d "$outdir/$name/IRMA/JSCDC_Segment"){
	`mkdir $outdir/$name/IRMA/JSCDC_Segment`;
}
for my $freq (@freq){
	my $freq_tmp = $freq * 100;
	my $freq_str = "consensusFreq".$freq_tmp;
	my $cons_dir = "$outdir/$name/IRMA/JSCDC_Segment/$freq_str";
	my $IRMA_module_dir = "HIV1freq".$freq_tmp."_JSCDC_Segment-pgm-hq"; # HIV1freq20_JSCDC_Segment
	print O "\# /data/fulongfei/git_repo/flu-amd/IRMA $IRMA_module_dir $fq $cons_dir\n\n";
}

# JSCDC_FullLength
if (!-d "$outdir/$name/IRMA/JSCDC_FullLength"){
	`mkdir $outdir/$name/IRMA/JSCDC_FullLength`;
}
for my $freq (@freq){
	my $freq_tmp = $freq * 100; # 2/5/10/15/20
	my $freq_str = "consensusFreq".$freq_tmp; # consensusFreq2/consensusFreq5/consensusFreq10/consensusFreq15/consensusFreq20
	my $cons_dir = "$outdir/$name/IRMA/JSCDC_FullLength/$freq_str";
	my $IRMA_module_dir = "HIV1freq".$freq_tmp."-pgm-hq"; # HIV1freq2/HIV1freq5/HIV1freq10/...
	print O "\# /data/fulongfei/git_repo/flu-amd/IRMA $IRMA_module_dir $fq $cons_dir\n\n";
}

print O "\# nextalign generateConsensus fasta\n";
my $cons_dir = "$outdir/$name/generateConsensus";
my $nextAlign_dir = "$outdir/$name/nextalign";
print O "\# perl /data/fulongfei/git_repo/HIVDrug/code/v1/merge_generateCons_fasta.pl $cons_dir $name $nextAlign_dir\n";
my $input_fa = "$outdir/$name/nextalign/input.fa";
print O "\# /data/fulongfei/git_repo/ClustalO_2019nCoV/bin/nextalign --sequences\=$input_fa --reference\=$ref --output-dir\=$outdir/$name/nextalign --output-basename\=$name --include-reference --in-order\n\n";
my $aln_fa = "$outdir/$name/nextalign/$name\.aligned.fasta";
my $parsed_file = "$outdir/$name/nextalign/$name\.parsed.txt";
print O "\# perl /data/fulongfei/git_repo/HIVDrug/scripts/parse_nextalign.pl $aln_fa $parsed_file\n\n";


print O "\# compare pileup and TVC freq\n";
my $pileup2vaf_variants_xls = "$outdir/$name/pileup2vaf/$name\.variants.xls";
my $cmp_txt = "$outdir/$name/freq_cmp.txt";
print O "\# perl /data/fulongfei/analysis/hiv/final_tvc_json_test/compare_variantCaller_and_pileup2vaf.pl $alleles_xls $pileup2vaf_variants_xls >$cmp_txt\n";

close O;

`chmod 755 $runsh`;





# http://10.69.138.39/report/275/#SARS_CoV_2_variantCaller-section

# Command : /results/plugins/SARS_CoV_2_variantCaller/bin/variant_caller_pipeline.py
	# --input-bam "/results/analysis/output/Home/Auto_user_GSS5PR-0043-179-Zheyi_416_275/plugin_out/SARS_CoV_2_variantCaller_out.215/IonCode_0109/IonCode_0109_rawlib.realigned.bam"
	# --region-bed "/results/uploads/BED/107/ion_ampliseq_sars-cov-2-insight/merged/plain/Ion_AmpliSeq_SARS-CoV-2-Insight.20210329.designed.bed"
	# --primer-trim-bed "/results/uploads/BED/107/ion_ampliseq_sars-cov-2-insight/unmerged/detail/Ion_AmpliSeq_SARS-CoV-2-Insight.20210329.designed.bed"
	# --generate-gvcf on
	# --postprocessed-bam "/results/analysis/output/Home/Auto_user_GSS5PR-0043-179-Zheyi_416_275/plugin_out/SARS_CoV_2_variantCaller_out.215/IonCode_0109/IonCode_0109_rawlib.realigned_processed.bam"
	# --reference-fasta "/results/referenceLibrary/tmap-f3/ion_ampliseq_sars-cov-2-insight/ion_ampliseq_sars-cov-2-insight.fasta"
	# --output-dir "/results/analysis/output/Home/Auto_user_GSS5PR-0043-179-Zheyi_416_275/plugin_out/SARS_CoV_2_variantCaller_out.215/IonCode_0109"
	# --parameters-file "/results/analysis/output/Home/Auto_user_GSS5PR-0043-179-Zheyi_416_275/plugin_out/SARS_CoV_2_variantCaller_out.215/local_parameters.json"
	# --bin-dir "/results/plugins/SARS_CoV_2_variantCaller/bin"
	# --error-motifs-dir "/results/plugins/SARS_CoV_2_variantCaller/share/TVC/sse"


# Command : /results/plugins/SARS_CoV_2_variantCaller/scripts/generate_variant_tables.py
	# --suppress-no-calls on
	# --input-vcf /results/analysis/output/Home/Auto_user_GSS5PR-0043-179-Zheyi_416_275/plugin_out/SARS_CoV_2_variantCaller_out.215/IonCode_0109/TSVC_variants.vcf
	# --region-bed "/results/uploads/BED/107/ion_ampliseq_sars-cov-2-insight/unmerged/detail/Ion_AmpliSeq_SARS-CoV-2-Insight.20210329.designed.bed"
	# --output-xls /results/analysis/output/Home/Auto_user_GSS5PR-0043-179-Zheyi_416_275/plugin_out/SARS_CoV_2_variantCaller_out.215/IonCode_0109/variants.xls
	# --alleles2-xls /results/analysis/output/Home/Auto_user_GSS5PR-0043-179-Zheyi_416_275/plugin_out/SARS_CoV_2_variantCaller_out.215/IonCode_0109/alleles.xls
	# --summary-json /results/analysis/output/Home/Auto_user_GSS5PR-0043-179-Zheyi_416_275/plugin_out/SARS_CoV_2_variantCaller_out.215/IonCode_0109/variant_summary.json
	# --scatter-png /results/analysis/output/Home/Auto_user_GSS5PR-0043-179-Zheyi_416_275/plugin_out/SARS_CoV_2_variantCaller_out.215/IonCode_0109/scatter.png
	# --barcode IonCode_0109
	# --concatenated-xls "/results/analysis/output/Home/Auto_user_GSS5PR-0043-179-Zheyi_416_275/plugin_out/SARS_CoV_2_variantCaller_out.215/R_2021_12_18_10_14_37_user_GSS5PR-0043-179-Zheyi.xls"
	# --run-name "R_2021_12_18_10_14_37_user_GSS5PR-0043-179-Zheyi"
	# --library-type "AmpliSeq"


# Processing IonCode_0109 (Sample 1)...

# (Sat Dec 18 13:56:52 CST 2021) Analyzing target contigs...
# Target regions: /results/uploads/BED/107/ion_ampliseq_sars-cov-2-insight/unmerged/detail/Ion_AmpliSeq_SARS-CoV-2-Insight.20210329.designed.bed
# References 1 contigs and 4 controls.
# (Sat Dec 18 13:56:52 CST 2021) Generating FASTA from VCF...
# /results/plugins/generateConsensus/gvcf_to_fasta.py
	# -m 1
	# -n 0.5
	# -p 0.6
	# -v '/results/analysis/output/Home/Auto_user_GSS5PR-0043-179-Zheyi_416_275/plugin_out/SARS_CoV_2_variantCaller_out.215/IonCode_0109/TSVC_variants.genome.vcf'
	# -o '/results/analysis/output/Home/Auto_user_GSS5PR-0043-179-Zheyi_416_275/plugin_out/generateConsensus_out.217/IonCode_0109/IonCode_0109_consensus.fasta'
	# -a 'Sample 1__2019-nCoV'
	# -c '2019-nCoV'
	# -d 20
	# -r /results/uploads/BED/107/ion_ampliseq_sars-cov-2-insight/unmerged/detail/Ion_AmpliSeq_SARS-CoV-2-Insight.20210329.designed.bed
