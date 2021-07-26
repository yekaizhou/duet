import logging
import argparse
import sys
import os

def set_logging(home):
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt = '%H:%M:%S')
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)
    fileHandler = logging.FileHandler(home + '/run_duet.log', mode="w")
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description = 'Duet: a fast and accurate germline variant caller for structural variation calling and haplotyping from low-depth nanopore sequencing reads.')
    parser.add_argument('-t', '--thread', type = int, default = 40,
        help = 'number of threads to use [%(default)s]')
    parser.add_argument('-m', '--min_allele_frequency', type = float, default = 0.25,
        help = 'minimum allele frequency required to call a candidate SNP [%(default)s]')
    parser.add_argument('-c', '--cluster_max_distance', type = float, default = 0.9,
        help = 'maximum span-position distance between SV marks in a cluster to call a SV candidates [%(default)s]')
    parser.add_argument('-s', '--sv_min_size', type = int, default = 50,
        help = 'minimum SV size to be reported [%(default)s]')
    parser.add_argument('-r', '--min_support_read', type = int, default = 2,
        help = 'minimum number of reads that support a SV to be reported [%(default)s]')
    parser.add_argument('BAM', type = str,
            help = 'sorted alignment file in .bam format (along with .bai file in the same directory)')
    parser.add_argument('REFERENCE', type = str,
            help = 'indexed reference genome in .fasta format (along with .fai file in the same directory)')
    parser.add_argument('OUTPUT', type = str,
        help = 'working and output directory (existing files in the directory will be overwritten)')
    args = parser.parse_args()
    return args

def check_envs(ref_path, aln_path):
    if not os.path.exists(aln_path + '.bai'):
        sys.exit("[ERROR] Alignment index .bai file not found, please run 'samtools index " + aln_path + "' first")
    if not os.path.exists(ref_path + '.fai'):
        sys.exit("[ERROR] Reference index .fai file not found, please run 'samtools faidx " + ref_path + "' first")