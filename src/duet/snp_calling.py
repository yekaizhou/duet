# coding=utf-8

import logging
import os
import time

def snp_calling(home, ref_path, aln_path, maf, thread, include_all_ctgs):
    lines = '*************************'
    logging.info(lines + ' SNP CALLING STARTED ' + lines)
    starttime = time.time()
    snp_calling_home = home + '/snp_calling/'
    os.system('mkdir ' + snp_calling_home)
    cmd = 'run_clair3.sh -b ' + aln_path + ' -f ' + ref_path + \
          ' -m "${CONDA_PREFIX}/bin/models/ont" -t ' + str(thread) + ' -p ont -o ' + snp_calling_home + \
          ' --snp_min_af=' + str(maf) + ' --pileup_only --call_snp_only'
    if include_all_ctgs:
        cmd = cmd + ' --include_all_ctgs'
    os.system(cmd)
    logging.info(lines + ' SNP CALLING COMPLETED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)