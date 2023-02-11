#!/bin/bash
  
workdir=~/
scripts_dir=$workdir/script
outdir=$workdir/out

R_dir=~/soft/R/v4.0.4/bin

ref_dir=~/working/reference
repeats=$ref_dir/ucsc/mm10/mm10_repeat_masker.gtf
transcriptome=$ref_dir/cellranger/refdata-gex-mm10-2020-A/genes/genes.gtf

cellranger_output=$workdir/data/cellranger_out

mkdir -p  $outdir

cd $workdir

# constructing spliced and unspliced counts matrices

#velocyto run10x -m $repeats $cellranger_output $transcriptome


# convert data from R/Seurat to Python/anndata
mkdir -p $outdir/data4anadata
$R_dir/Rscript $scripts_dir/cellrank0_data4anadata.R    $outdir 1>$outdir/cellrank0.log

python3 $scripts_dir/cellrank1_seurat2anndata.py 1>$outdir/cellrank1.log 2>&1

# 
python3 $scripts_dir/cellrank2_merge_Seurat_and_velocyto.py 1>$outdir/cellrank2.log 2>&1

#
mkdir -p $outdir/RNA_velocity
python3 $scripts_dir/cellrank3_RNA_velocity_analysis.py 1>$outdir/cellrank3.log 2>&1

