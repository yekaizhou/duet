import logging
import time
from duet.write_file import print_sv_header, print_sv
from duet.sv_phasing_fn import generate_phased_callset

def sv_phasing(home, svlen_thres, suppread_thres, thread):
    lines = '*************************'
    logging.info(lines + ' SV PHASING STARTED ' + lines)
    starttime = time.time()
    sv_calling_path = home + '/sv_calling/variants.vcf'
    sv_phasing_path = home + '/phased_sv.vcf'
    snp_phasing_home = home + '/snp_phasing/'
    logging.info('create output .vcf file')
    print_sv_header(sv_calling_path, sv_phasing_path)
    print_sv(generate_phased_callset(sv_calling_path, snp_phasing_home, svlen_thres, suppread_thres, thread), \
             sv_phasing_path)
    logging.info(lines + ' SV PHASING COMPLETED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)