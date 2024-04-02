#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p cn
#SBATCH -J Marina_SComatic
#SBATCH -o SComatic_Pancreas_Marina.o
source $HOME/.bashrc
bash_dir=/public/home/lorihan/Shell/shell.database
db=/public/home/lorihan/lrh/hg38/refdata-gex-GRCh38-2020-A
db_pre=/public/home/lorihan/lrh/hg38/GRCh38-2020_premrna
ls $db
fq_dir=/public/home/lorihan/lrh/database/Pancreas/Marina/FastQ
# mkdir $fq_dir/../cellranger_count

# cat $fq_dir/../Marina_pdata.txt|while read id
# do
# name=`echo $id | awk -F ';' '{print $1}'`
# cellranger_outDir=`echo $id | awk -F ';' '{print $22}'`

# fq1=$fq_dir/${name}_1.fastq.gz
# fq2=$fq_dir/${name}_2.fastq.gz

# if [ -e $fq2 ]; then
# echo $cellranger_outDir $name 
# ls $fq2 -htlr
# mv -i  $fq1 $fq_dir/${cellranger_outDir}_S1_L001_R1_001.fastq.gz
# mv -i  $fq2 $fq_dir/${cellranger_outDir}_S1_L001_R2_001.fastq.gz

# cd $fq_dir/../cellranger_count
# echo "start Single Cell cellranger for $cellranger_outDir" `date`
# if [ ! -e $cellranger_outDir ]; then
# echo $cellranger_outDir $new
# cellranger count --id=$cellranger_outDir \
# --localcores=32 \
# --localmem=64 \
# --transcriptome=$db \
# --fastqs=$fq_dir \
# --sample=$cellranger_outDir \
# --chemistry auto
# # --expect-cells=9000 --nosecondary(maybe timesaving)
# fi
# fi
# done

# Prepare files and environment
# Activate conda environment if needed
source $HOME/.bashrc
conda activate SComatic
# Create an output folder and go to the main SComatic folder
data_path=/public/home/lorihan/lrh/database/Pancreas/Marina
SCOMATIC=/public/home/lorihan/lrh/Evolution/SCI-supp/SComatic
output_dir=$data_path/SCOMATIC
mkdir -p $output_dir
outdir=$data_path/scMapping
mkdir $outdir
# If reference genome is not unpacked and indexed (.fai), you have to do it before running SComatic
cd $output_dir
REF=/public/home/lorihan/lrh/hg38/refdata-gex-GRCh38-2020-A/fasta/genome.fa

# Step 1: Splitting alignment file in cell type specific bams
cat $data_path/Marina_pdata.txt|while read id
do
name=`echo $id | awk -F ';' '{print $1}'`
sample=`echo $id | awk -F ';' '{print $22}'`
if [ ! -e $sample ]; then
mkdir $sample
echo $name $sample
echo "start call SNV for $sample" `date`
# sample='NAG_1'
output_dir1=$output_dir/$sample/Step1_BamCellTypes
mkdir -p $output_dir1

# mv $output_dir/../cellranger_count/${sample}/outs/${sample}*possorted_genome_bam.bam.bai  $output_dir/../cellranger_count/${sample}/outs/possorted_genome_bam.bam.bai
# samtools index $output_dir/../cellranger_count/${sample}/outs/possorted_genome_bam.bam && \
# echo -e "\n*** possorted_genome_bam.bam index done at $(date +'%T %F') ***\n" 

python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py --bam $output_dir/../cellranger_count/${sample}/outs/possorted_genome_bam.bam \
        --meta $data_path/cell_barcode_annotations.tsv \
        --id ${sample} \
        --n_trim 5 \
        --max_nM 5 \
        --max_NH 1 \
        --outdir $output_dir1
# Step 2: Collecting base count information

output_dir2=$output_dir/$sample/Step2_BaseCellCounts
mkdir -p $output_dir2

for bam in $(ls -d $output_dir1/*bam);do
  
  # Cell type
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')

  # Temp folder
  temp=$output_dir2/temp_${cell_type}
  mkdir -p $temp

  # Command line to submit to cluster
  python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
    --ref $REF \
    --chrom all \
    --out_folder $output_dir2 \
    --min_bq 30 \
    --tmp_dir $temp \
    --nprocs 32

  rm -rf $temp
done
# In our experience, when running SComatic in an HPC cluster it is most efficient to compute the base count information for all cell types at once, specially when submitting each python line independently and at the same time to the cluster, and using multiple processors (p.e. --nprocs 32)

# Step 3: Merging base count matrices
# sample=Example
output_dir3=$output_dir/$sample/Step3_BaseCellCountsMerged
mkdir -p $output_dir3

python $SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder ${output_dir2} \
  --outfile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv
# Step 4: Detection of somatic mutations
# Step 4.1 & Step 4.2
# Step 4.1
output_dir4=$output_dir/$sample/Step4_VariantCalling
mkdir -p $output_dir4

# sample=Example
# REF=$SCOMATIC/example_data/chr10.fa

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
          --infile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv \
          --outfile ${output_dir4}/${sample} \
          --ref $REF

# Step 4.2
editing=$SCOMATIC/RNAediting/AllEditingSites.hg38.txt
PON=$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py \
          --infile ${output_dir4}/${sample}.calling.step1.tsv \
          --outfile ${output_dir4}/${sample} \
          --editing $editing \
          --pon $PON
# As done in the datasets analysed in our study, we suggest intersecting the variants with the high quality regions of the human genome (provided in SComatic/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed). Hence, if you want to do this intersection and keep the somatic mutations that passed all filters, you can do it with a simple command:

bedtools intersect -header -a ${output_dir4}/${sample}.calling.step2.tsv -b $SCOMATIC/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed | awk '$1 ~ /^#/ || $6 == "PASS"' > ${output_dir4}/${sample}.calling.step2.pass.tsv
# Other SComatic functionalities
# Computing the number of callable sites per cell type
# sample=Example
output_dir5=$output_dir/$sample/CellTypeCallableSites
mkdir -p $output_dir5

python $SCOMATIC/scripts/GetCallableSites/GetAllCallableSites.py --infile $output_dir4/${sample}.calling.step1.tsv  \
   --outfile $output_dir5/${sample} \
   --max_cov 150 --min_cell_types 2
# Computing the number of callable sites per cell
# sample=Example
# REF=$SCOMATIC/example_data/chr10.fa
STEP4_1=$output_dir4/${sample}.calling.step1.tsv

output_dir6=$output_dir/$sample/UniqueCellCallableSites
mkdir -p $output_dir6

for bam in $(ls -d $output_dir1/*bam);do  
    cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
    echo $cell_type
    
    temp=$output_dir6/temp_${cell_type}
    mkdir -p $temp

    python $SCOMATIC/scripts/SitesPerCell/SitesPerCell.py --bam $bam    \
       --infile $output_dir4/${sample}.calling.step1.tsv   \
       --ref $REF \
       --out_folder $output_dir6 --tmp_dir $temp --nprocs 32
    echo
done
# Computing the genotype for each cell at the variant sites
META=$data_path/cell_barcode_annotations.tsv
# sample=Example
# REF=$SCOMATIC/example_data/chr10.fa
STEP4_2_pass=${output_dir4}/${sample}.calling.step2.pass.tsv

output_dir7=$output_dir/$sample/SingleCellAlleles
mkdir -p $output_dir7

for bam in $(ls -d $output_dir1/*bam);do  
    cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
    
    temp=$output_dir7/temp_${cell_type}
    mkdir -p $temp

    python $SCOMATIC/scripts/SingleCellGenotype/SingleCellGenotype.py --bam $bam  \
        --infile ${STEP4_2_pass}   \
        --nprocs 32   \
        --meta $META   \
        --outfile ${output_dir7}/${cell_type}.single_cell_genotype.tsv  \
        --tmp_dir $temp  \
        --ref $REF

    rm -rf $temp
done

echo "start annoFilter for $sample " `date`
cd $outdir
output_dir8=$outdir/$sample
mkdir -p $output_dir8
META=$data_path/cell_barcode_annotations.tsv
cat ${output_dir4}/${sample}.calling.step2.tsv|grep 'PASS' > $output_dir8/${sample}.PASS.tsv
STEP4_2_pass=$output_dir8/${sample}.PASS.tsv

for bam in $(ls -d $output_dir1/*bam);do  
cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')

temp=$output_dir8/temp_${cell_type}
mkdir -p $temp

python $SCOMATIC/scripts/SingleCellGenotype/SingleCellGenotype.py --bam $bam  \
--infile ${STEP4_2_pass}   \
--nprocs 32   \
--meta $META   \
--outfile ${output_dir8}/${cell_type}.single_cell_genotype.tsv  \
--tmp_dir $temp  \
--ref $REF

rm -rf $temp
done

fi
done
