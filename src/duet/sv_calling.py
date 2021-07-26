import logging
import os
import time

def sv_calling(home, ref_path, aln_path, cls_thres, svlen_thres):
    lines = '*************************'
    logging.info(lines + ' SV CALLING STARTED ' + lines)
    starttime = time.time()
    sv_calling_home = home + '/sv_calling/'
    os.system('mkdir ' + sv_calling_home)
    os.system('svim alignment ' + sv_calling_home + ' ' + aln_path + ' ' + ref_path + ' --min_sv_size ' + \
              str(svlen_thres) + ' --read_names --minimum_depth 0 --minimum_score 0 --cluster_max_distance ' + \
              str(cls_thres))
    logging.info(lines + ' SV CALLING COMPLETED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)