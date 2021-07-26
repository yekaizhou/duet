# coding=utf-8

import logging
import os
import time

def snp_calling(home, ref_path, aln_path, maf, thread):
    lines = '*************************'
    logging.info(lines + ' SNP CALLING STARTED ' + lines)
    starttime = time.time()
    snp_calling_home = home + '/snp_calling/'
    os.system('mkdir ' + snp_calling_home)
    os.system('run_clair3.sh -b ' + aln_path + ' -f ' + ref_path + \
              ' -m "${CONDA_PREFIX}/bin/models/ont" -t ' + str(thread) + ' -p ont -o ' + snp_calling_home + \
              ' --snp_min_af=' + str(maf) + ' --pileup_only --call_snp_only')
    logging.info(lines + ' SNP CALLING COMPLETED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)