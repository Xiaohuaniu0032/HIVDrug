ref='/data/fulongfei/git_repo/HIVDrug/IRMA/K03455.fasta'
bed='/data/fulongfei/git_repo/HIVDrug/IRMA/JSCDC.Pol.Gene.bed'
bedtools getfasta -fi $ref -bed $bed -fo consensus.fasta
