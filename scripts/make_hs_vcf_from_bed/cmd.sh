bed='/data/fulongfei/analysis/hiv/re_analysis_new_ref/POL.bed'
bed='HIV.drug.HS.BED'
fa='/data/fulongfei/git_repo/HIVDrug/ref/K03455.fasta'

/data/fulongfei/tools/variantCaller/bin/tvcutils prepare_hotspots --input-bed $bed --reference $fa --left-alignment on  --allow-block-substitutions on --output-bed "$PWD/POL.left.bed" --output-vcf "$PWD/POL.hotspot.vcf" --unmerged-bed /data/fulongfei/git_repo/HIVDrug/BED/new_BED/WG00471_HIV_pol.20201027_Designed.NewRef.bed 
