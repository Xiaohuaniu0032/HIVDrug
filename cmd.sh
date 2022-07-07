#perl get_drug_DNA_and_AA.pl $PWD >drug.aa.list.xls
gvcf='/data/fulongfei/analysis/hiv/re_analysis_new_ref_with_hs/IonXpress_002_rawlib/variantCaller/TSVC_variants.genome.vcf'
name='IonXpress_002_rawlib'
ref='/data/fulongfei/git_repo/HIVDrug/ref/K03455.fasta'
bed='/data/fulongfei/analysis/hiv/re_analysis_new_ref/POL.bed'
bed='/data/fulongfei/git_repo/HIVDrug/BED/new_BED/WG00471_HIV_pol.20201027_Designed.NewRef.bed'
outdir="$PWD/test"

perl HIVDrug.pl $gvcf $name $ref $bed $outdir
