import subprocess
import time

def main():
    #sizes = [0.4, 0.8, 1.6, 3.2, 6.4]
    # only need to run for 3.2 and 6.4, already ran for 0.4, 0.8, 1.6
    sizes = [3.2, 6.4]
    num_runs = 30

    for size in sizes:
        # generate all metagenome filenames
        metagenome_files = [ "metagenomes/metagenome_"+str(size)+"_seed_"+str(seed)+".fastq" for seed in range(1, num_runs+1) ]
        
        # run diamond batch for all these metagenems
        # command: python ../../run_tools/run_diamond/run_diamond_batch.py --threads 128 --diamond_script ../../run_tools/run_diamond/run_diamond.py --output_dir diamond_fast_batch_output <metagenomes>
        # also, record time and memory

        start_time = time.time()

        cmd = "/usr/bin/time -v python ../../run_tools/run_diamond/run_diamond_batch.py --threads 128 --diamond_script ../../run_tools/run_diamond/run_diamond.py --output_dir diamond_fast_batch_output " + " ".join(metagenome_files)
        print(cmd)
        output = subprocess.check_output(cmd, stderr = subprocess.STDOUT, shell=True)
        output = output.decode("utf-8")
        output = output.split('\n')

        end_time = time.time()

        print(output)

        # user time is at line 23 from the bottom
        # system time is at line 22 from the bottom
        # maximum resident set size is at line 15 from the bottom
        # average resident set size is at line 14 from the bottom
        # cpu time is the sum of user time and system time
        
        user_time = float(output[-23].split(" ")[-1].strip())
        system_time = float(output[-22].split(" ")[-1].strip())
        cpu_time = user_time + system_time
        max_resident_set_size = float(output[-15].split(" ")[-1].strip())
        avg_resident_set_size = float(output[-14].split(" ")[-1].strip())
        walltime = end_time - start_time

        # write to file: single row: walltime,cputime,maxresident,avgresident
        with open("diamond_fast_batch_output/diamond_resource_usage_" + str(size) + ".txt", "w") as f:
            f.write(str(walltime) + "," + str(cpu_time) + "," + str(max_resident_set_size) + "," + str(avg_resident_set_size))

if __name__ == '__main__':
    main()