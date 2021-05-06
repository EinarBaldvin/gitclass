#!/bin/bash

set -euxo pipefail

## genotype list - one ID per line

list=samples.to.extract_reduced
genes=genes_cov_95perc.exonum_lead0.bed
vcf_snp=var_snps/extracted-16.snpeff_gff.all.478sam.dp8_gq20.HC.recal99.9.newtruth.wo_fs.SNP.biall.var_nonvar1.chr.vcf
vcf_nonvar=novar_snps/extract.vcf
vcf_indel=indel_snps/478sam.dp8_gq20.HC.filtered.INDEL.biall.var_nonvar1.chr.vcf
#horvu=$1
#exon=$2
#chr=$3
#start=$4
#end=$5
#strand=$6

# Read in list of genes to be extracted
# v5 change start and end of extraction
# I hope never to use this again!!!!
while read -r chr start end strand horvu exon

do

if [ -z "$list" ]; then #1

	cat $vcf_snp | awk '$1 == '${chr}' && $2 >= '${start}' && $2 <= '${end}'' > ${horvu}.${exon}.snp

	cat $vcf_indel | awk '$1 == '${chr}' && $2 >= '${start}' && $2 <= '${end}'' > ${horvu}.${exon}.indel

	cat $vcf_nonvar | awk '$1 == '${chr}' && $2 >= '${start}' && $2 <= '${end}'' > ${horvu}.${exon}.nonvar

	grep CHR -m 1 $vcf_snp | cat - <(cat ${horvu}.${exon}.snp ${horvu}.${exon}.indel ${horvu}.${exon}.nonvar | sort -k2,2n) > temp.${horvu}.${exon}.all


	vcf_fin=temp.${horvu}.${exon}.all


else

cut_snp=`cat ${list} | cut -f1 | while read line; do awk '{for(i=1;i<=NF;i++)if($i == "'${line}'")a=i}{print a}' <(grep "#CHR" -m 1 ${vcf_snp}); done | tr -s "\n" "," | sed 's:,$::g' | awk '{print"cut -f1-9,"$0" '${vcf_snp}'"}'`

cat <(eval ${cut_snp} | grep "CHR" -m 1) <(eval ${cut_snp} | awk '$1 == '${chr}' && $2 >= '${start}' && $2 <= '${end}'') > temp.var.${list}

cut_indel=`cat ${list} | cut -f1 | while read line; do awk '{for(i=1;i<=NF;i++)if($i == "'${line}'")a=i}{print a}' <(grep "#CHR" -m 1 ${vcf_indel}); done | tr -s "\n" "," | sed 's:,$::g' | awk '{print"cut -f1-9,"$0" '${vcf_indel}'"}'`

eval ${cut_indel} | awk '$1 == '${chr}' && $2 >= '${start}' && $2 <= '${end}'' > temp.indel.${list}

cut_nonvar=`cat ${list} | cut -f1 | while read line; do awk '{for(i=1;i<=NF;i++)if($i == "'${line}'")a=i}{print a}' <(grep "#CHR" -m 1 ${vcf_nonvar}); done | tr -s "\n" "," | sed 's:,$::g' | awk '{print"cut -f1-9,"$0" '${vcf_nonvar}'"}'`

eval ${cut_nonvar} | awk '$1 == '${chr}' && $2 >= '${start}' && $2 <= '${end}'' > temp.nonvar.${list}

cat temp.nonvar.${list} temp.indel.${list} temp.var.${list} | sort -u -k2,2n > temp.${horvu}.${exon}.all.${list}

vcf_temp=temp.${horvu}.${exon}.all.${list}



fi #1

# Trying to fix missing refs

for i in {10..22}; do cut -f$i $vcf_temp | sed 's/:PASS//g' | awk 'BEGIN{FS=":|"; OFS=":"}{if(NF <= 4 && $1 ~ "\\./\\." && $3 == 0 && $2 > 30)sub("\\./\\.","0/0",$1);print}' | sed 's/:0:/,0:/' > $i.temp; done

paste <(cut -f 1-9 $vcf_temp) {10..22}.temp > temp.fixed.${horvu}.${exon}.all.${list}

rm {10..22}.temp

vcf_fin=temp.fixed.${horvu}.${exon}.all.${list}

wait

num=`grep "#CHR" $vcf_fin -m1  | awk '{print NF}'`


for ((i=10;i<=${num};i++)); do 
	id=`grep "#CHR" $vcf_fin -m1  | 
		cut -f${i}`;  
		grep -v "#" ${vcf_fin} | 
		cut -f4,5,${i} |  
		awk '{if(length($1)+length($2) == 2 && $3 ~ "^0/0") print substr($1,1,1); 
	else if(length($1)+length($2) == 2 && $3 ~ "^1/1") print substr($2,1,1);  
	else if(length($2) > 1 && length($1) == 1 &&  $3 ~ "^0/0") print substr($1,1,1);
	else if(length($2) > 1 && length($1) == 1 &&  $3 ~ "^1/1") {gsub(/,.*/,"",$2); print $2; printf ("\n");}
	else if(length($2) == 1 && length($1) > 1 &&  $3 ~ "^0/0") print substr($1,1,1);
	else if(length($2) == 1 && length($1) > 1 &&  $3 ~ "^1/1") print substr($2,1,1); 
	else if(length($2) > 1 && length($1) > 1 && $3 ~ "^0/0") print substr($1,1,1);
	else if(length($2) > 1 && length($1) > 1 && $3 ~ "^1/1") {gsub(/,.*/,"",$2); print $2; printf ("\n");}
	else print ""}' | 
		tr -d "\n" | 
		sed '$s/$/\n/' | 
		cat <(echo ">"$id" "${horvu}" "${exon}) -;
	done > temp.${horvu}.${exon}.fa 

if [ $strand == "-" ]; then

revseq -auto -stdout temp.${horvu}.${exon}.fa | fasta_formatter -w 0 | sed 's/Reversed: /rev_comp_/' > ${horvu}_${exon}.rev.fa

rm temp.${horvu}.${exon}.fa

else

mv temp.${horvu}.${exon}.fa ${horvu}_${exon}.fa

fi

#rm $vcf_fin

done < ${genes}

wait

# Paste everything together, copy ref CDS fasta from Morex !!! taken out because of using MACSEv2 (stop codons removed at the end if present)!!! This is such a mess!
cat ${genes} | awk '{print $5}' | sort -u | while read line; do cat <(paste -d "\0" $line\_*) <(grep -m 1 -A 1 $line ../reference/160517_Hv_IBSC_PGSB_r1_CDS_HighLowConf_REPR_annotation.unwrapped.fasta | sed 's/...$//' ) | sed '/ 01.*//' > ext_genes/$line.fa ; done

wait

# Remove all temp files
cat ${genes} | awk '{print $5}' | sort -u | while read line; do rm "$line"_* ; done
