#!/usr/bin/env python3
# coding=utf-8

import os
import logging
import time
import sys
import numpy as np
from operator import itemgetter
from duet.utils import set_logging, parse_args, check_envs
# from duet.snp_calling import snp_calling
from duet.snp_phasing import snp_phasing
from duet.sv_calling import sv_calling
# from duet.sv_phasing import sv_phasing
from duet.sv_phasing_fn import generate_callinfo, read_hap_bam
from duet.read_file import init_chrom_list
from duet.write_file import *

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

def get_phase_info(call, ps_sr = 2, ps_num = 1, oneps_set = ''):
    info = dict()
    info['hap1'] = info['hap2'] = info['tothap'] = info['hap1_totsc'] = info['allhap'] \
    = info['hap2_totsc'] = info['hap1_avgsc'] = info['hap2_avgsc'] = info['ps'] = info['hap0'] = 0
    if ps_num == 1:
        for read in call['svreadinfo']:
            if len(read) > 1:
                info['ps'] = read[2]
                if read[1] == 1:
                    info['hap1'] += 1
                    info['hap1_totsc'] += read[3]
                elif read[1] == 2:
                    info['hap2'] += 1
                    info['hap2_totsc'] += read[3]
        info['hap0'], info['allhap'] = 0, info['hap1'] + info['hap2']
    if ps_num == 2:
        ps_dict = dict()
        for read in call['svreadinfo']:
            if len(read) > 1:
                info['allhap'] += 1
                ps = read[2]
                if ps in oneps_set:
                    if ps not in ps_dict:
                        ps_dict[ps] = dict()
                        ps_dict[ps][0] = ps_dict[ps][1] = ps_dict[ps][2] = ps_dict[ps]['1_score'] = \
                        ps_dict[ps]['2_score'] = 0
                    ps_dict[ps][read[1]] += 1
                    ps_dict[ps][str(read[1]) + '_score'] += read[3]
                    ps_dict[ps][0] += 1
        max_tot = 0
        for ps in ps_dict:
            if ps_dict[ps][0] > max_tot:
                max_tot = ps_dict[ps][0]
                info['hap1'], info['hap1_totsc'], info['hap2'], info['hap2_totsc'], info['ps'] = \
                ps_dict[ps][1], ps_dict[ps]['1_score'], ps_dict[ps][2], ps_dict[ps]['2_score'], ps
                info['hap0'] = info['allhap'] - info['hap1'] - info['hap2']
    if ps_num == 0 or info['hap1'] == info['hap2'] == 0:
        oneps_ls = np.sort(np.array(list(oneps_set)))
        idx = np.searchsorted(oneps_ls, call['pos']) # map pos to ls
        idx = max(idx - 1, 0) if abs(call['pos'] - oneps_ls[max(idx - 1, 0)]) < \
              abs(call['pos'] - oneps_ls[min(idx, len(oneps_ls) - 1)]) else min(idx, len(oneps_ls) - 1)
        info['ps'] = oneps_ls[idx]
    info['hapread_ratio'] = info['allhap'] / len(call['svreadinfo'])
    if info['hap1'] > 0:
        info['hap1_avgsc'] = info['hap1_totsc'] / info['hap1']
    if info['hap2'] > 0:
        info['hap2_avgsc'] = info['hap2_totsc'] / info['hap2']
    if info['hap1'] >= ps_sr:
        info['tothap'] += 1
    if info['hap2'] >= ps_sr:
        info['tothap'] += 2
    info['nohap'] = len(call['svreadinfo']) - info['allhap']
    info['hap_diff'] = abs(info['hap1'] - info['hap2'])
    info['sv_ratio'] = call['svread'] / (call['svread'] + call['refread'])
    info['totsc_ratio'] = max(info['hap1_totsc'], info['hap2_totsc']) / min(info['hap1_totsc'], info['hap2_totsc']) \
                          if min(info['hap1_totsc'], info['hap2_totsc']) > 0 else 0
    info['onehap_totsc'] = max(info['hap1_totsc'], info['hap2_totsc']) if \
                           min(info['hap1_totsc'], info['hap2_totsc']) == 0 else 0
    info['avgsc_ratio'] = max(info['hap1_avgsc'], info['hap2_avgsc']) / min(info['hap1_avgsc'], info['hap2_avgsc']) \
                          if min(info['hap1_avgsc'], info['hap2_avgsc']) > 0 else 0
    info['onehap_avgsc'] = max(info['hap1_avgsc'], info['hap2_avgsc']) if \
                           min(info['hap1_avgsc'], info['hap2_avgsc']) == 0 else 0
    info['hap_avgsc_diff'] = abs(info['hap2_avgsc'] - info['hap1_avgsc'])
    info['hap_totsc_diff'] = abs(info['hap2_totsc'] - info['hap1_totsc'])
    info['ref_num'] = call['refread']
    info['sv_num'] = call['svread']
    info['allsv'] = len(call['svreadinfo'])
    info['hap_ratio'] = max(info['hap1'], info['hap2']) / max(min(info['hap1'], info['hap2']), 1)
    info['totsc'] = info['hap1_totsc'] + info['hap2_totsc']
    info['totsc_ratio2'] = max(info['hap1_totsc'], info['hap2_totsc']) / info['totsc'] if info['totsc'] > 0 else 0
    return info

def generate_phased_callset(vcf_path, sam_home, svlen_thres, suppread_thres, thread):
    callstat = generate_callinfo(vcf_path, read_hap_bam(sam_home, thread))
    logging.info('integrate read weight information')
    chrom_list = init_chrom_list()
    callset = [c for c in callstat if c['svlen'] >= svlen_thres and \
               c['svread'] >= suppread_thres and c['callgt'] not in ['./.']]
    callset_dict = dict()
    callset_dict[1] = [c for c in callset if len(set([s[2] for s in c['svreadinfo'] if len(s) > 1])) == 1]
    callset_dict[2] = [c for c in callset if len(set([s[2] for s in c['svreadinfo'] if len(s) > 1])) > 1]
    callset_dict[0] = [c for c in callset if len(set([s[2] for s in c['svreadinfo'] if len(s) > 1])) == 0]
    oneps_set = [set() for ch in range(24)]
    logging.info('calculate read weight statistics')
    for ch in range(24):
        callset_ch = [c for c in callset_dict[1] if c['chrom'] in ['chr' + chrom_list[ch], chrom_list[ch]]]
        for c in callset_ch:
            for read in c['svreadinfo']:
                if len(read) > 1:
                    oneps_set[ch].add(read[2])
                    break
    phased_callset = []
    logging.info('predict SV haplotypes in the callset')
    for ch in range(24):
        for ps_num in range(1, 2):
            callset_ch = [c for c in callset_dict[ps_num] if c['chrom'] in ['chr' + chrom_list[ch], chrom_list[ch]]]
            for c in callset_ch:
                pred, ps = predict_hp(c, ps_num, oneps_set[ch])
                if pred == 0:
                    continue
                phased_callset.append(dict())
                phased_callset[-1]['ps'] = ps
                if pred == 1:
                    phased_callset[-1]['hp'] = '1|0'
                if pred == 2:
                    phased_callset[-1]['hp'] = '0|1'
                if pred == 3:
                    phased_callset[-1]['hp'] = '1|1'
                phased_callset[-1]['chrom'] = c['chrom']
                phased_callset[-1]['pos'] = c['pos']
                phased_callset[-1]['svlen'] = c['svlen'] if c['svtype'] == 'INS' else -c['svlen']
                phased_callset[-1]['alt'] = c['svtype']
    phased_callset.sort(key = itemgetter('chrom', 'pos'))
    return phased_callset

def predict_hp(call, ps_num, oneps_set):
    f = get_phase_info(call, 2, ps_num, oneps_set)
    pred = 0
    if ps_num != 1:
        return pred, f['ps']
    else:
        return f['tothap'], f['ps']

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