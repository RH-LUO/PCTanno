#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p cn
source $HOME/.bashrc
human_index_dir=/data/rluo4/lorihan/hg38/refdata-gex-GRCh38-2020-A
# mouse_index_dir=/data/rluo4/lorihan/hg38/refdata-gex-mm10-2020-A
fastqs_dir=/data/rluo4/database/Cervix/Guo/space_ranger
image_path=/data/rluo4/database/Cervix/Guo/Anterior_image.jpg
output_dir=/data/rluo4/database/Cervix/Guo/spaceranger_count
mkdir $output_dir; mkdir $image_path
bash_dir=/data/rluo4/lorihan/Shell/

ls $human_index_dir
cd $output_dir
cat /data/rluo4/All/Output/GEO/Guo/Guo_ST_pdata.txt|sort -nr|while read id
do
name=`echo $id | awk -F ';' '{print $4}'`
sample=`echo $id | awk -F ';' '{print $2}'`
fq1=$fastqs_dir/${name}_1.fastq.gz
fq2=$fastqs_dir/${name}_2.fastq.gz
if [ ! -e $sample ]; then
# if [ -e $fq2 ]; then
echo $sample $name
echo $fq1
# mv -i  $fq1 $fastqs_dir/${sample}_S1_L001_R1_001.fastq.gz
# mv -i  $fq2 $fastqs_dir/${sample}_S1_L001_R2_001.fastq.gz
#>>>A.sh>>>
cd ${output_dir}
#Anterior样本定量
spaceranger count \
    --id $sample \
    --description Cervix_Anterior_Section_1 \
    --transcriptome ${human_index_dir} \
    --fastqs ${fastqs_dir} \
    --image ${image_path}/${sample}*.jpg \
    --slide V19L29-035 \
    --area B1 \
    --sample $sample \
    --localcores 20 \
    --localmem 128
#<<<A.sh<<<
# fi
fi
done
