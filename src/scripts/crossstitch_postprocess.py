import shlex
import subprocess
import numpy as np
import argparse

def init_chrom_list():
    chrom_list = []
    for i in range(1, 23):
        chrom_list.append(str(i))
    chrom_list.append('X')
    chrom_list.append('Y')
    return(chrom_list)

thread = 40

parser = argparse.ArgumentParser(description = 'extract SV calls from crossstitch for benchmarking')
parser.add_argument('cs_dir', type = str, help = 'crossstitch working and output directory')
crossstitch_home = parser.parse_args().cs_dir

vcf_header = subprocess.check_output(shlex.split('bcftools view -h --threads ' + str(thread) + ' ' + crossstitch_home + '/cs.spliced.scrubbed.vcf')).decode('ascii').split('\n')[:-1]
crossstitch_vcf_str = subprocess.check_output(shlex.split('bcftools view -H --threads ' + str(thread) + ' ' + crossstitch_home + '/cs.spliced.scrubbed.vcf')).decode('ascii').split('\n')[:-1]
crossstitch_vcf_str = [s.split() for s in crossstitch_vcf_str]
chrom_list = init_chrom_list()
ps_list = [[] for ch in range(24)]
for ch in range(24):
    vcf_ch = [s for s in crossstitch_vcf_str if s[0] == 'chr' + chrom_list[ch]]
    snp_ch = [s for s in vcf_ch if s[2] == '.']
    if not snp_ch:
        continue
    if snp_ch[0][-1].split(':')[3].find('/') > 0:
        ps_list[ch] = sorted([int(s) for s in list(set([s[-1].split(':')[2] for s in snp_ch]) - {'.'})])
    else:
        ps_list[ch] = sorted([int(s) for s in list(set([s[-1].split(':')[4] for s in snp_ch if len(s[-1].split(':')) == 5]) - {'.'})])
tot_vcf_str = ''
for ch in range(24):
    vcf_ch = [s for s in crossstitch_vcf_str if s[0] == 'chr' + chrom_list[ch]]
    sv_ch = [s for s in vcf_ch if s[2] != '.']
    phased_sv_ch = [s for s in sv_ch if '|' in s[-1]]
    for sv in phased_sv_ch:
        idx = np.searchsorted(ps_list[ch], int(sv[1]))
        if idx == len(ps_list[ch]) or (idx > 0 and abs(ps_list[ch][idx] - int(sv[1])) > abs(ps_list[ch][idx - 1] - int(sv[1]))):
            idx -= 1 
        tot_vcf_str += '\t'.join(sv[:-2]) + '\tHP:PS\t' + sv[-1][:3] + ':' + str(ps_list[ch][idx]) + '\n'
dip_vcf = open(crossstitch_home + '/phased_sv.vcf', 'w')
for s in vcf_header:
    dip_vcf.write(s + '\n')
dip_vcf.write(tot_vcf_str)
dip_vcf.close()
