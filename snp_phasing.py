import logging
import os
import time
import subprocess
import shlex

def snp_phasing(home, ref_path, aln_path, thread):
    lines = '*************************'
    logging.info(lines + ' SNP PHASING STARTED ' + lines)
    starttime = time.time()
    snp_calling_home = home + '/snp_calling/'
    snp_phasing_home = home + '/snp_phasing/'
    os.system('mkdir ' + snp_phasing_home)
    pileup_vcf_path = snp_calling_home + 'pileup.vcf.gz'
    ctgs = subprocess.check_output(shlex.split('tabix --list-chroms ' + pileup_vcf_path)).decode('ascii').split('\n')[:-1]
    for chrs in ctgs:
        os.system('bcftools view -r ' + chrs + ' -c1 ' + pileup_vcf_path + ' > ' + snp_phasing_home + chrs + '.vcf')
    f = open(snp_phasing_home + 'parallel_wh.sh', 'w')
    f.write('CHR=(' + ' '.join(ctgs) + ')\n')
    f.write('parallel -j' + str(thread) + ' "whatshap phase -o ' + snp_phasing_home + 'phased_{1}.vcf.gz -r ' + \
            ref_path + ' --chromosome {1} --distrust-genotypes --ignore-read-groups ' + \
            snp_phasing_home + '{1}.vcf ' + aln_path + '" ::: ${CHR[@]}\n')
    f.write('parallel -j' + str(thread) + ' "tabix -f -p vcf ' + snp_phasing_home + \
            'phased_{1}.vcf.gz" ::: ${CHR[@]}\n')
    f.write('parallel -j' + str(thread) + ' "whatshap haplotag -o ' + snp_phasing_home + '{1}.bam -r ' + \
            ref_path + ' --regions {1} --ignore-read-groups --tag-supplementary ' + \
            snp_phasing_home + 'phased_{1}.vcf.gz ' + aln_path + '" ::: ${CHR[@]}\n')
    f.close()
    os.system('chmod a+x ' + snp_phasing_home + 'parallel_wh.sh')
    os.system(snp_phasing_home + 'parallel_wh.sh')
    logging.info(lines + ' SNP PHASING COMPLETED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)