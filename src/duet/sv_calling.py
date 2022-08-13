# coding=utf-8

import logging
import os
import time

def sv_calling(home, ref_path, aln_path, cls_thres, svlen_thres, thread, caller, supp_thres):
    lines = '*************************'
    logging.info(lines + ' SV CALLING STARTED ' + lines)
    starttime = time.time()
    sv_calling_home = home + '/sv_calling/'
    os.system('mkdir ' + sv_calling_home)
    if caller == 'svim':
        os.system('svim alignment ' + sv_calling_home + ' ' + aln_path + ' ' + ref_path + ' --min_sv_size ' + \
                  str(svlen_thres) + ' --read_names --minimum_depth 0 --minimum_score 0 --cluster_max_distance ' + str(cls_thres))
    else:
        os.system('cuteSV --genotype --report_readid ' + aln_path + ' ' + ref_path + ' ' + sv_calling_home + 'variants.vcf ' \
                  + sv_calling_home + ' -t ' + str(thread) + ' -s ' + str(supp_thres) + ' -l ' + str(svlen_thres))
    logging.info(lines + ' SV CALLING COMPLETED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)