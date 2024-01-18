import argparse
import datetime
import configparser
import subprocess

# hardcoded items
# designed for the GPU machine
#camisim_path = '/home/grads/mbr5797/camisim/CAMISIM-1.3/metagenome_from_profile.py'
camisim_path = '/home/grads/mbr5797/camisim/CAMISIM-1.3/metagenomesimulation.py'

def main():
    # parse arguments
    parser = argparse.ArgumentParser(description='Process number of genomes.')
    parser.add_argument('number_of_genomes', type=int, help='Number of genomes')
    parser.add_argument('--config', type=str, help='Config file')
    parser.add_argument('--seed', type=int, help='Seed used in simulation', default=0)
    parser.add_argument('--outdir', type=str, help='Output directory', default='./out')
    parser.add_argument('--size', type=float, help='Size of the file in Gbp', default=0.1)
    args = parser.parse_args()

    # read the arguments
    number_of_genomes = args.number_of_genomes
    if args.config:
        config_filename = args.config
    else:
        timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        config_filename = 'config_' + timestamp + '.ini'
    seed = args.seed
    outdir = args.outdir
    size = args.size

    # read the sample config file
    config = configparser.ConfigParser()
    config.read('sample_config.ini')

    # update the num_genomes field with the command line argument
    # section: Main
    config.set('Main', 'seed', str(seed))
    config.set('Main', 'output_directory', outdir)

    # section: ReadSimulator
    config.set('ReadSimulator', 'size', str(size))

    # section: community0
    config.set('community0', 'genomes_total', str(number_of_genomes))
    config.set('community0', 'genomes_real', str(number_of_genomes))

    # write the updated config to the config filename
    with open(config_filename, 'w') as config_file:
        config.write(config_file)

    # run camisim
    cmd = f'python {camisim_path} -seed {seed} -s 2 -id KL -p 64 {config_filename}'
    print(cmd)
    subprocess.run(cmd.split(' '))

if __name__ == '__main__':
    main()
