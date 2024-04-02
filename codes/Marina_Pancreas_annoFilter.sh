#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cn
#SBATCH -J Anno_Marina.nr
#SBATCH -o Anno_Pancreas_Marina.nr.o
source $HOME/miniconda3/bin/activate
data_path=/public/home/lorihan/lrh/database/Pancreas/Marina
SCOMATIC=/public/home/lorihan/lrh/Evolution/SCI-supp/SComatic
output_dir=$data_path/SCOMATIC
indir=$data_path/anno.var
mkdir $indir
GATK_bundle="$HOME/lrh/hg38"
ref=/public/home/lorihan/lrh/hg38/refdata-gex-GRCh38-2020-A/fasta/genome.fa
tool_path="$HOME/lrh/bin"
bin=$HOME/lrh/bin/annovar
db=$HOME/lrh/bin/annovar/humandb/
cat $data_path/Marina_pdata.txt|sort -nr|while read id
do
name=`echo $id | awk -F ';' '{print $1}'`
sample=`echo $id | awk -F ';' '{print $22}'`
if [ -e $output_dir/$sample ]; then
if [ ! -e $indir/${sample}.hg38_multianno.txt ]; then
output_dir4=$output_dir/$sample/Step4_VariantCalling
line=`cat $output_dir4/*.calling.step2.tsv|grep -v ^#|wc -l`
echo $sample $line 
cat $output_dir4/*.calling.step2.tsv | grep -v '^#'| grep 'PASS' |  tr '\t' '-' | awk -F'-' -v OFS='\t' '{print $1,$2,$3,$4,$5,$0}' >$indir/${sample}.variants.avinput
echo "start annoFilter for $sample " `date`
perl $bin/table_annovar.pl $indir/${sample}.variants.avinput \
                $db -buildver hg38 \
                -out $indir/${sample} \
                -remove -protocol refGene,cytoBand,genomicSuperDups,rmsk,avsnp150,dbnsfp42c,exac03,1000g2015aug_all,gnomad_genome -operation g,r,r,r,f,f,f,f,f \
                -nastring NA
fi
fi
done
