import os
import logging
import time
import sys
from utils import set_logging, parse_args, check_envs
from snp_calling import snp_calling
from snp_phasing import snp_phasing
from sv_calling import sv_calling
from sv_phasing import sv_phasing

def main(argv):
    args = parse_args(argv)
    check_envs(args.REFERENCE, args.BAM)
    if not os.path.exists(args.OUTPUT):
        os.makedirs(args.OUTPUT)
    set_logging(args.OUTPUT)
    lines = '*************************'
    starttime = time.time()
    logging.info(lines + ' DUET STARTED ' + lines)
    snp_calling(args.OUTPUT, args.REFERENCE, args.BAM, args.min_allele_frequency, args.thread)
    sv_calling(args.OUTPUT, args.REFERENCE, args.BAM, args.cluster_max_distance, args.sv_min_size)
    snp_phasing(args.OUTPUT, args.REFERENCE, args.BAM, args.thread)
    sv_phasing(args.OUTPUT, args.sv_min_size, args.min_support_read, args.thread)
    logging.info(lines + ' DUET FINISHED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)
    logging.info('OUTPUT .VCF FILE AT ' + args.OUTPUT + '/phased_sv.vcf')

if __name__ == '__main__':
    main(sys.argv[1:])