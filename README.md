# run_camisim
In this repo, we will write a python wrapper, that will run camisim simulations. The script will 

1. take all necessary arguments,
1. create a config file as required by camisim, and
1. invoke the camisim program.

# Details

This script runs camisim with the given number of genomes, 
    and then converts the bam files to ground truth files. Finally,
    it merges the fastq files into a single file.
    Inputs: 
    - number of genomes to simulate (int)
    - config file (will be written by this script, and used by camisim) (str)
    - seed for simulations (int)
    - output directory for camisim (str)
    - size of the file in Gbp (float)
    - gene ground truth filename (str)
    - ko ground truth filename (str)
    - merged metagenome filename (str)
    
# Usage example
```
python run_camisim.py 10 --config config_seed_0_size_0.1.ini --seed 0 --outdir ./out_seed_0_size_0.1 --size 0.1 --gene_g_t gene_ground_truth_seed_0_size_0.1.csv --ko_g_t ko_ground_truth_seed_0_size_0.1.csv --metagenome_filename metagenome_seed_0_size_0.1.fastq
```

# Next to do

1. Add code to run the tools: diamond, sourmash gather, and mmseqs2. The code should also generate their resource usages
1. Produce outputs of these tools
1. Write code to read the outputs and assess the performance
1. Write code to aggregate the resource usages

# Running simulations plan

Run varying

1. number of genomes [in progress]
1. size of the sample
1. different error profiles
1. number of novel genomes
1. divergence rate
1. other taxonomy
1. other parameters... (distribution)