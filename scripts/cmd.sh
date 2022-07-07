perl MakeHotspotVCF.pl ../drug.aa.list.xls /data/fulongfei/git_repo/HIVDrug/db $PWD
perl MakeHotspotBed.pl HIV.drug.HS.vcf >HIV.drug.HS.BED

perl parse_nextalign.pl ../test/IonXpress_002_rawlib.aligned.fasta >align.o
