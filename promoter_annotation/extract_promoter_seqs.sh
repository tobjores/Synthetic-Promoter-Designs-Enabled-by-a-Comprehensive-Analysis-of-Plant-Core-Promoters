#!/bin/bash

# requires:
#	- bedtools (version 2.29.2)
#	- bedops (version 2.4.35)
#	- R (version 4.0.0) with libraries:
#		- dplyr (version 1.0.0)
#		- tidyr (version 1.1.0)
#		- readr (version 1.3.1)
#		- tibble (version 3.0.1)


### download genomes and annotations ###
if [ ! -f "Arabidopsis_genome.fa" ]; then
  curl https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas > Arabidopsis_genome.fa
fi
if [ ! -f "Arabidopsis_annotation.gff3" ]; then
  curl https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz | gunzip | grep -v "^#" | sed 's/^Chr//g' > Arabidopsis_annotation.gff3
fi

if [ ! -f "Maize_genome.fa" ]; then
  curl ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz | gunzip > Maize_genome.fa
fi
if [ ! -f "Maize_annotation.gff3" ]; then
  curl ftp://ftp.ensemblgenomes.org/pub/plants/release-42/gff3/zea_mays/Zea_mays.B73_RefGen_v4.42.chr.gff3.gz | gunzip | grep -v "^#" > Maize_annotation.gff3
fi

if [ ! -f "Sorghum_genome.fa" ]; then
  curl ftp://ftp.ensemblgenomes.org/pub/plants/release-43/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz | gunzip > Sorghum_genome.fa
fi
if [ ! -f "Sorghum_annotation.gff3" ]; then
  curl ftp://ftp.ensemblgenomes.org/pub/plants/release-43/gff3/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.43.chr.gff3.gz | gunzip | grep -v "^#" > Sorghum_annotation.gff3
fi


### get promoter sequences for Arabidopsis, maize, and sorghum ###
for SPECIES in "Arabidopsis" "Maize" "Sorghum"
do

	### extract coordinates of the chromosomal genes ###
	if [ "${SPECIES}" = "Arabidopsis" ]
	then
		awk '{if($3 == "gene" && $1 ~ /[0-9]+/) print}' ${SPECIES}_annotation.gff3 \
			| grep "type=mirna" \
			| awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' \
			> ${SPECIES}_miRNA_genes.gff3
	elif [ "${SPECIES}" = "Maize" ] || [ "${SPECIES}" = "Sorghum" ]
	then
		awk '{if($3 == "ncRNA_gene" && $1 ~ /[0-9]+/) print}' ${SPECIES}_annotation.gff3 \
			| grep "type=pre_miRNA" \
			| awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' \
			> ${SPECIES}_miRNA_genes.gff3
	fi

	awk '{if($3 == "gene" && $1 ~ /[0-9]+/) print}' ${SPECIES}_annotation.gff3 \
		| grep "type=protein_coding" \
		| awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' \
		> ${SPECIES}_protein_coding_genes.gff3


	### generate bed file with promoter coordinates (+5 to -165) ###
	if [ ${SPECIES} = 'Maize' ]
	then
		Rscript correct_Maize_TSSs.R
	fi
	
	awk -v OFS='\t' '{
		ID=substr($9, 4, index($9, ";") - 4);
		sub(/gene:/, "", ID);
		if($7 == "+") {print $1, ($4 - 166), ($4 + 4), ID, 0, $7}
		else if($7 == "-") {print $1, ($5 - 5), ($5 + 165), ID, 0, $7}
	}' ${SPECIES}_protein_coding_genes.gff3 | awk '$2 >= 0' > ${SPECIES}_protein_coding_promoters.bed

	awk -v OFS='\t' '{
		ID=substr($9, 4, index($9, ";") - 4);
		sub(/gene:/, "", ID);
		if($7 == "+") {print $1, ($4 - 166), ($4 + 4), ID, 0, $7}
		else if($7 == "-") {print $1, ($5 - 5), ($5 + 165), ID, 0, $7}
	}' ${SPECIES}_miRNA_genes.gff3 | awk '$2 >= 0' > ${SPECIES}_miRNA_promoters.bed


	### get fasta sequences ###
	bedtools getfasta -s -nameOnly -fi ${SPECIES}_genome.fa -bed ${SPECIES}_protein_coding_promoters.bed > ${SPECIES}_protein_coding_promoters.fasta
	
	bedtools getfasta -s -nameOnly -fi ${SPECIES}_genome.fa -bed ${SPECIES}_miRNA_promoters.bed > ${SPECIES}_miRNA_promoters.fasta

	
	### mutate BsaI (GGTCTC/GAGACC)/BbsI (GAAGAC/GTCTTC) sites and list mutations (T>A in BsaI fwd; A>T in BsaI rev; G>C in BbsI fwd; C>G in BbsI rev) ###
	awk -v OFS='\t' 'BEGIN{print "gene", "type", "sequence", "strand", "mutations"} {
	  GENE=substr($1, 2, index($1, "(") - 2); STRAND=substr($1, index($1, "(") + 1, 1);
	  getline;
	  if($1 !~ /N/) {
		SEQ=$1; MUT="";
		while(index(SEQ, "GGTCTC") > 0) {MUT=MUT index(SEQ, "GGTCTC") + 2 "T>A;"; sub(/GGTCTC/, "GGACTC", SEQ)};
		while(index(SEQ, "GAGACC") > 0) {MUT=MUT index(SEQ, "GAGACC") + 3 "A>T;"; sub(/GAGACC/, "GAGTCC", SEQ)};
		while(index(SEQ, "GAAGAC") > 0) {MUT=MUT index(SEQ, "GAAGAC") + 3 "G>C;"; sub(/GAAGAC/, "GAACAC", SEQ)};
		while(index(SEQ, "GTCTTC") > 0) {MUT=MUT index(SEQ, "GTCTTC") + 2 "C>G;"; sub(/GTCTTC/, "GTGTTC", SEQ)};
		sub(/;$/, "", MUT);
		print GENE, "protein_coding", SEQ, STRAND, MUT
	  }
	}' ${SPECIES}_protein_coding_promoters.fasta > ${SPECIES}_all_promoters.tsv

	awk -v OFS='\t' '{
	  GENE=substr($1, 2, index($1, "(") - 2); STRAND=substr($1, index($1, "(") + 1, 1);
	  getline;
	  if($1 !~ /N/) {
		SEQ=$1; MUT="";
		while(index(SEQ, "GGTCTC") > 0) {MUT=MUT index(SEQ, "GGTCTC") + 2 "T>A;"; sub(/GGTCTC/, "GGACTC", SEQ)};
		while(index(SEQ, "GAGACC") > 0) {MUT=MUT index(SEQ, "GAGACC") + 3 "A>T;"; sub(/GAGACC/, "GAGTCC", SEQ)};
		while(index(SEQ, "GAAGAC") > 0) {MUT=MUT index(SEQ, "GAAGAC") + 3 "G>C;"; sub(/GAAGAC/, "GAACAC", SEQ)};
		while(index(SEQ, "GTCTTC") > 0) {MUT=MUT index(SEQ, "GTCTTC") + 2 "C>G;"; sub(/GTCTTC/, "GTGTTC", SEQ)};
		sub(/;$/, "", MUT);    
		print GENE, "miRNA", SEQ, STRAND, MUT
	  }
	}' ${SPECIES}_miRNA_promoters.fasta >> ${SPECIES}_all_promoters.tsv

	
	### mutate BsaI/BbsI sites forming with the 5' or 3' cloning adapter ###
	awk -v OFS='\t' '{
	  if ($3 ~ /^GTCTC/) {SEQ=$3; sub(/^GTCTC/, "GACTC", SEQ); MUT="2T>A;" $5; print $1, $2, SEQ, $4, MUT}
	  else if ($3 ~ /^AGACC/) {SEQ=$3; sub(/^AGACC/, "AGTCC", SEQ); MUT="3A>T;" $5; print $1, $2, SEQ, $4, MUT}
	  else if ($3 ~ /^AAGAC/) {SEQ=$3; sub(/^AAGAC/, "AACAC", SEQ); MUT="3G>C;" $5; print $1, $2, SEQ, $4, MUT}
	  else if ($3 ~ /^TCTTC/) {SEQ=$3; sub(/^TCTTC/, "TGTTC", SEQ); MUT="2C>G;" $5; print $1, $2, SEQ, $4, MUT}
	  else print $0
	}' ${SPECIES}_all_promoters.tsv > tmp

	mv tmp ${SPECIES}_all_promoters.tsv

	if [ ${SPECIES} = 'Arabidopsis' ]
	then
	  awk -v OFS='\t' '{
		if ($3 ~ /GG$/) {SEQ=$3; sub(/GG$/, "GA", SEQ); MUT=$5 length(SEQ) "G>A"; print $1, $2, SEQ, $4, MUT}
		else if ($3 ~ /GGTC$/) {SEQ=$3; sub(/GGTC$/, "GGAC", SEQ); MUT=$5 length(SEQ) - 1 "T>A"; print $1, $2, SEQ, $4, MUT}
		else if ($3 ~ /GTCT$/) {SEQ=$3; sub(/GTCT$/, "GTGT", SEQ); MUT=$5 length(SEQ) - 1 "C>G"; print $1, $2, SEQ, $4, MUT}
		else print $0
	  }' ${SPECIES}_all_promoters.tsv > tmp

	  mv tmp ${SPECIES}_all_promoters.tsv
	elif [ ${SPECIES} = 'Maize' ]
	then
	  awk -v OFS='\t' '{
		if ($3 ~ /GGTCT$/) {SEQ=$3; sub(/GGTCT$/, "GGACT", SEQ); MUT=$5 length(SEQ) - 2 "T>A"; print $1, $2, SEQ, $4, MUT}
		else if ($3 ~ /GAGAC$/) {SEQ=$3; sub(/GAGAC$/, "GAGTC", SEQ); MUT=$5 length(SEQ) - 1 "A>T"; print $1, $2, SEQ, $4, MUT}
		else if ($3 ~ /GAAGA$/) {SEQ=$3; sub(/GAAGA$/, "GAACA", SEQ); MUT=$5 length(SEQ) - 1 "G>C"; print $1, $2, SEQ, $4, MUT}
		else if ($3 ~ /GTCTT$/) {SEQ=$3; sub(/GTCTT$/, "GTGTT", SEQ); MUT=$5 length(SEQ) - 2 "C>G"; print $1, $2, SEQ, $4, MUT}
		else print $0
	  }' ${SPECIES}_all_promoters.tsv > tmp

	  mv tmp ${SPECIES}_all_promoters.tsv
	fi


	### identify protein-coding genes without an annotated 5' UTR (TSS = first base of CDS) ###
	awk -v FS='[\t;]' -v OFS='\t' '$3 == "CDS" {sub("ID=", "", $9); print $1, $4 - 1, $5, $9, $6, $7}' ${SPECIES}_annotation.gff3 \
		| sort-bed - > ${SPECIES}_CDS.bed

	awk -v FS='[\t;]' -v OFS='\t' '$4 < $5 {sub("ID=(gene:)?", "", $9); print $1, $4 - 1, $5, $9, $6, $7}' ${SPECIES}_protein_coding_genes.gff3 \
		| sort-bed - > ${SPECIES}_protein_coding_genes.bed

	bedmap --echo --echo-map-range --delim '\t' ${SPECIES}_protein_coding_genes.bed ${SPECIES}_CDS.bed \
		| awk -v OFS='\t' -v species=${SPECIES} '{
			print $1, $8, $9, $4, $5, $6;
			if (($6 == "+" && $2 == $8) || ($6 == "-" && $3 == $9)) print $4 > species"_noUTR.txt"
		}' > ${SPECIES}_CDS-range.bed
	
done


### Correct a "K" nucleotide in the promoter of AT2G39050 ###
sed -i 's/K/T/g' Arabidopsis_all_promoters.tsv


### collapse identical sequences and create reference sequences ###
Rscript unique_promoters.R