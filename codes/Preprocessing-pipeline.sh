#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p cn
# 1) cellranger
source $HOME/.bashrc
bash_dir=$HOME/Shell/shell.database
db=$HOME/lrh/hg38/refdata-gex-mm10-2020-A
ls $db
fq_dir=$HOME/database/Pancreas/Zhang/FastQ
mkdir $fq_dir/../cellranger_count
cd $fq_dir/../cellranger_count
cat $fq_dir/../Zhang_pdata.txt|while read id
do
name=`echo $id | awk -F ';' '{print $4}'`
sample=`echo $id | awk -F ';' '{print $1}'`
fq1=$fq_dir/${name}_1.fastq.gz
fq2=$fq_dir/${name}_2.fastq.gz
fq3=$fq_dir/${name}_3.fastq.gz

cellranger count --id=$sample \
--localcores=16 \
--localmem=64 \
--transcriptome=$db \
--fastqs=$fq_dir \
--sample=$sample \
--chemistry auto

done

# 2) velocyto
source $HOME/miniconda3/bin/activate
conda activate pyvelo
data_path="$HOME/lrh/database/Pancreas/Zhang"
GATK_bundle="$HOME/lrh/hg38"
rmsk_gtf="$GATK_bundle/mm10_rmsk.gtf"
cellranger_output=${data_path}/cellranger_count
cellranger_gtf=$HOME/lrh/hg38/refdata-gex-mm10-2020-A/genes/genes.gtf
ls -lh $rmsk_gtf  $cellranger_output $cellranger_gtf

cd $cellranger_output
ls  */outs/possorted_genome_bam.bam|while read id
do  new=${id/possorted_genome_bam.bam/cellsorted_possorted_genome_bam.bam}
echo $new  $id
if [ ! -e $new ]; then
echo "start samtools for $id " `date`
samtools sort -@ 16  -t CB -O BAM -o $new  $id
fi
done

cat $data_path/Zhang_pdata.txt|while read id
do
new=cellsorted_possorted_genome_bam.bam
name=`echo $id | awk -F ';' '{print $4}'`
cellranger_outDir=`echo $id | awk -F ';' '{print $1}'`

cd $cellranger_output
if [ -e $cellranger_outDir/outs/$new ]; then
echo $name $cellranger_outDir
# mv -i $cellranger_outDir/velocyto $cellranger_outDir/velocyto_incorrect
if [ ! -e $cellranger_outDir/velocyto ]; then
mv -i $cellranger_outDir/outs/filtered_feature_bc_matrix $cellranger_outDir/outs/filtered_feature_bc_matrix_my
mv $data_path/Ori_barcode/$cellranger_outDir $cellranger_outDir/outs/filtered_feature_bc_matrix
cd $data_path
echo "start velocyto for $cellranger_outDir" `date`
velocyto run10x $cellranger_output/$cellranger_outDir $cellranger_gtf -m $rmsk_gtf
fi
fi
done