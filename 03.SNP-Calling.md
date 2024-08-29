# 02. Stacks Pipeline

## SNAKEFILE

This [**02.snakefile_stacks**]() is designed to process BAM files using the Stacks software suite. It includes steps for running the ref_map pipeline, populations analysis, correcting allelic imbalance, and performing PCA and MDS analyses on the resulting data. The output files are generated in a structured directory according to your configuration.

You shouldn’t need to modify anything in the snakefile unless you want to tweak parameters that are not explicitly set in the config file already! If you want to run the pipeline as is, just modify the [**02.snakefile_stacks.yml**]() config file with the appropriate file paths.

![Stacks Pipeline DAG](02.Stacks-DAG.png)

## REQUIRED ACCESSORY SCRIPTS

Part of this pipeline also relies on the following accessory scripts:

* [**allelicBalance.py**]() - For correcting allelic imbalance in the VCF file.
* [**append_population.py**]() - For appending population information to the VCF file.
* [**pca_plot.R**]() - For visualizing PCA results.
* [**mds_plot.R**]() - For visualizing MDS results.

You'll want to make sure you download these scripts and keep them in the same directory as your Snakemake and .yml files.

## OTHER NECESSARY INPUT FILES

* **Merged BAM Files**: The pipeline expects BAM files to be located in the directory specified by the `merge_dir` path in your config file.
* **Population Map**: A file called [**popmap.txt**]() containing two columns (without headers). The first column is the sample ID, and the second column is the population ID.
* **Indexed Reference Genome**: Ensure that your reference genome is indexed and the file path is correctly set in the configuration.
* **Reference Nickname**: A short nickname for your reference genome to simplify file naming and management.

## PIPELINE CONFIGURATION

The configuration file [**02.snakefile_stacks.yml**]() should be updated with paths and parameters specific to your project:

- **Project Directory**: Update the `home_dir` path to point to your project's root directory.
- **Populations Run Name**: Specify the name of your populations run in `populations_run`.
- **Reference Genome**: Ensure the `reference_file` and `reference_nickname` are correctly set.
- **Population Parameters**: Adjust `min_maf`, `max_obs_het`, and `min_samples_overall` as needed.
- **Allelic Balance Correction**: Set the `p_value` and `ratio` for allelic balance correction.
- **BCFTools Parameters**: Modify the `missingness` parameter if necessary.

## RUNNING THE PIPELINE

Once you have all the necessary files and have modified the file paths (and any settings you want to change) in the [**02.snakefile_stacks.yml**]() config file, you can run the pipeline with the following commands:

If you haven't already created a conda environment and configured your SLURM profile, see [**00.Conda-Snakemake-Slurm**](/00.Conda-Snakemake-Slurm).

```bash
# I like to run this using screen
screen -S StacksWorkflow

# Load the conda environment
micromamba activate GBS

# This will do a quick dry run of your pipeline.
snakemake -s 02.snakefile_stacks --profile slurm -n -r

# If everything looks good (lots of green and yellow printout — no red), you can execute the pipeline. Adjust the number of jobs as needed.
snakemake -s 02.snakefile_stacks --profile slurm -j 20
