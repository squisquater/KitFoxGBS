#!/bin/bash -l
#SBATCH --job-name=CanisVulpesUrocyonAncestralFasta
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 3:00:00
#SBATCH --mem=4GB
#SBATCH -p high
#SBATCH -o /group/ctbrowngrp2/sophiepq/KitFoxGBS/angsd/SFS/slurmlogs/CanisVulpesUrocyonAncestralFasta.out
#SBATCH -e /group/ctbrowngrp2/sophiepq/KitFoxGBS/angsd/SFS/slurmlogs/CanisVulpesUrocyonAncestralFasta.err

module load samtools/1.10
module load angsd

cd /group/ctbrowngrp2/sophiepq/KitFoxGBS/angsd/SFS/AncestralState
POP_BAMLIST=/group/ctbrowngrp2/sophiepq/KitFoxGBS/angsd/SFS/AncestralState/CanisVulpesUrocyonAncestralBamlist.txt
OUTNAME=CanisVulpesUrocyonAnc

angsd -bam $POP_BAMLIST -dofasta 2 -P 10 -out $OUTNAME -setMinDepth 60 -doCounts 1 -setMaxDepth 1000 -minMapQ 1 -minQ 20 -remove_bads 1 -uniqueOnly 1
