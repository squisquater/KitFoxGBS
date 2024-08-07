# 01. Data Processing

The following are steps to convert your raw sequencing data to sample specific FASTQ files that can then be aligned to a reference genome for downstream analysis. 

*Note: The PE merged fastq files resulting from this Data Processing Pipeline can be found in the NCBI SRA under BioProject No. PRJXXXXXXXXX. If you need the raw data please reach out to squisquater@ucdavis.edu.*

## Step 1: Demultiplex Data

**Description**:
**`process_radtags`** is designed to process raw FASTQ data from RAD-tag sequencing protocols, like GBS. It demultiplexes the input FASTQ files based on barcode sequences and optionally cleans the data to remove low-quality reads.

- The barcode file needs to be in the following format (no headers, column1 = unique barcode, column2 = sample ID

```python
AACT	S12-0649
CCAG	S12-0650
AACCA	S12-0651
CCACG	S12-0653
TATAA	S12-0654
GAGCG	S12-0655
ACATA	S12-0656
...   ...
...   ...
```

**Parameters**:

- **`-1 <path_to_R1>`**: Path to the first (forward) paired-end FASTQ file.
- **`-2 <path_to_R2>`**: Path to the second (reverse) paired-end FASTQ file.
- **`-o <output_directory>`**: Path to the directory where demultiplexed output files should be written.
- **`-b <barcode_file>`**: Path to the file containing barcode sequences used for demultiplexing.
- **`-c`**: Clean the data to remove reads with uncalled bases.
- **`-q`**: Remove low-quality reads. The program will discard reads with an average quality below a threshold (usually Phred score 10, but can be adjusted with additional parameters).
- **`-r`**: Rescue barcodes and RAD-tags. This will correct barcodes and RAD-tags that are within one Levenshtein distance from a true barcode (allows for 1 mismatch).
- **`-e <enzyme>`**: Specify the restriction enzyme used during RAD or GBS library preparation.
    
    Example: **`ecoT22I`**
    

**Example Shell Script**

*update all filenames/filepaths accordingly*

```python
nano KFGBS1_processradtags_PE.sh
```

```bash
#!/bin/bash -l
#SBATCH --job-name=KFGBS1_processradtags_PE
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH -t 24:00:00
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH --mem=5G
#SBATCH -o "/group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914/slurmlogs/KFGBS1_processradtags_PE.out"
#SBATCH -e "/group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914/slurmlogs/slurmlogs/KFGBS1_processradtags_PE.err"

STARTTIME=$(date +"%s")

micromamba activate GBS #you'll need to swap out micromamba for mamba or conda depending on what you installed.

process_radtags -1 /group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914/Undetermined_S0_L004_R1_001.fastq.gz -2 /group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914/Undetermined_S0_L004_R2_001.fastq.gz -o /group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914/Demultiplexed_PE -b /group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914/KitFox_GBS1_barcodefile.txt -c -q -r -e ecoT22I

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Took $timestamp hours:minutes:seconds to complete..."
```

```python
sbatch KFGBS1_processradtags_PE.sh # this submits the job!
```

## Step 2: Trim Data

**Description**:

- This script uses the **`trim_galore`** tool to:
    - Trim low-quality bases and adapter sequences from the paired-end FASTQ files.
    - Run FastQC quality checks on the trimmed sequences.
- To run this script ensure you have a file named **`samples.txt`** in your data directory (**`DIR`**) with each line corresponding to a sample ID. The script will process samples based on the order of this file in conjunction with the **`SLURM_ARRAY_TASK_ID`**. The script is currently set up to run 96 samples with a max of 32 samples at a time( **`--array=1-96%32`**).
- Set the **`DIR`** variable to the path where the raw GBS data (.fastq files) is located.

**Example Shell Script**

*update all filenames/filepaths accordingly*

```python
nano KFGBS1_trimPE.sh
```

```bash
#!/bin/bash -l
#SBATCH --job-name=KFGBS1_trimPE
#SBATCH --array=1-96%32
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH -t 45:00
#SBATCH -p bmh
#SBATCH -A ctbrowngrp
#SBATCH --mem=500M
#SBATCH -o /group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914/slurmlogs/KFGBS1_trimPE_%A_%a.out
#SBATCH -e /group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914/slurmlogs/KFGBS1_trimPE_%A_%a.err

STARTTIME=$(date +"%s")

micromamba activate GBS #you'll need to swap out micromamba for mamba or conda depending on what you installed.

DIR="/group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914"
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${DIR}/samples.txt | cut -f1)

trim_galore --fastqc --paired ${DIR}/Demultiplexed_PE/${SAMPLE}.1.fq.gz ${DIR}/Demultiplexed_PE/${SAMPLE}.2.fq.gz --output_dir ${DIR}/trimmed

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Took $timestamp hours:minutes:seconds to complete..."
```

```python
sbatch KFGBS1_trimPE.sh
```

## Step 3: Merge Data

**Description**:

- This script utilizes the PEAR (Paired-End reAd mergeR) software to merge paired-end reads from previously trimmed FASTQ files. After merging, it organizes the resultant files by moving them to a designated directory.
- Ensure that you have previously trimmed FASTQ files present in the **`${DIR}/trimmed`** directory with the naming pattern **`${SAMPLE}.1_val_1.fq.gz`** and **`${SAMPLE}.2_val_2.fq.gz`**.

**Example Shell Script**

*update all filenames/filepaths accordingly*

```python
nano KFGBS1_mergePE.sh
```

```bash
#!/bin/bash -l
#SBATCH --job-name=KFGBS1_mergePE
#SBATCH --array=1-96%32
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH -t 45:00
#SBATCH -p bmh
#SBATCH -A ctbrowngrp
#SBATCH --mem=500M
#SBATCH -o /group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914/slurmlogs/KFGBS1_mergePE_%A_%a.out
#SBATCH -e /group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914/slurmlogs/KFGBS1_mergePE_%A_%a.err

STARTTIME=$(date +"%s")

micromamba activate GBS #you'll need to swap out micromamba for mamba or conda depending on what you installed.

DIR="/group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914"
##shouldn't need to modify anything below this line
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${DIR}/samples.txt | cut -f1)

pear -f ${DIR}/trimmed/${SAMPLE}.1_val_1.fq.gz -r ${DIR}/trimmed/${SAMPLE}.2_val_2.fq.gz -o ${SAMPLE}

mv ${SAMPLE}.assembled.fastq ${DIR}/merged
mv ${SAMPLE}.discarded.fastq ${DIR}/merged
mv ${SAMPLE}.unassembled.forward.fastq ${DIR}/merged
mv ${SAMPLE}.unassembled.reverse.fastq ${DIR}/merged

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Took $timestamp hours:minutes:seconds to complete..."
```

```python
sbatch KFGBS1_mergePE.sh
```

## Step 4: Run FastQC

**Description**:

- This command runs the FastQC tool on the merged FASTQ file of a specified sample. FastQC provides a quality check on raw sequence data, producing various charts and metrics that can be helpful in understanding the data's quality.
- You can download the resulting .html files and view them in your web browser.

```bash
#!/bin/bash -l
#SBATCH --job-name=KFGBS1_fastqcPE.sh
#SBATCH --array=1-96%32
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH -t 45:00
#SBATCH -p bmh
#SBATCH -A ctbrowngrp
#SBATCH --mem=500M
#SBATCH -o /group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914/slurmlogs/KFGBS1_fastqcPE_%A_%a.out
#SBATCH -e /group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914/slurmlogs/KFGBS1_fastqcPE_%A_%a.err

STARTTIME=$(date +"%s")

micromamba activate GBS #you'll need to swap out micromamba for mamba or conda depending on what you installed.

DIR="/group/ctbrowngrp2/sophiepq/RawGBSData/KFGBS1_20220914"
##shouldn't need to modify anything below this line
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${DIR}/samples.txt | cut -f1)

fastqc ${DIR}/merged/${SAMPLE}.assembled.fastq

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Took $timestamp hours:minutes:seconds to complete..."
```

