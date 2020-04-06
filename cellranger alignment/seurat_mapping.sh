#!/bin/bash -l

id=$1
shift
data=$1
shift
sample=$1
shift

#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:15:00     # 1 day and 15 minutes
#SBATCH --output=$out.fastqc.stdout
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="cellranger"
#SBATCH -p highmem # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

#Load the packages
/bigdata/messaoudilab/sloanal/Projects/Alcohol_BM/scRNAseq_CD34/CellRanger/cellranger-3.1.0/cellranger count --id $id --transcriptome=/bigdata/messaoudilab/sloanal/Projects/Alcohol_BM/scRNAseq_CD34/ref/Mmul_8 --fastq=$data --sample=$sample --localcores=32 --expect-cells=12000
