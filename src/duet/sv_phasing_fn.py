# coding=utf-8

import logging
import os
import shlex
import subprocess
import numpy as np
from operator import itemgetter
from read_file import init_chrom_list, parse_vcf

def read_hap_bam(path, thread, include_all_ctgs):
    logging.info('extract SNP signatures')
    chrom_list = init_chrom_list(include_all_ctgs, path[:len(path) - 13])
    
    read_hap = [[] for ch in range(len(chrom_list))]
    
    for ch in range(len(chrom_list)):
        read_hap[ch] = dict()
        if os.path.exists(path + 'chr' + chrom_list[ch] + '.bam'):
            hap_bam_path = path + 'chr' + chrom_list[ch] + '.bam'
        elif os.path.exists(path + chrom_list[ch] + '.bam'):
            hap_bam_path = path + chrom_list[ch] + '.bam'
        else:
            continue
        alns = subprocess.check_output(shlex.split('samtools view -@' + str(thread) + ' ' + hap_bam_path)).decode('ascii').split('\n')[:-1]
        for s in alns:
            s = s.split()
            if 'PC:i:' in s[-2]:
                read_hap[ch][s[0]] = {'hap': int(s[-3][5:]), 'ps': int(s[-1][5:]), 'pc': int(s[-2][5:])}
        if alns:
            logging.info('  signatures extracted from ' + chrom_list[ch])
        else:
            logging.info('  no signature from ' + chrom_list[ch])
    return read_hap

def generate_callinfo(caller_path, read_hap, include_all_ctgs):
    logging.info('extract SV signatures')
    comp_call = parse_vcf(caller_path, include_all_ctgs)
    chrom_list = init_chrom_list(include_all_ctgs, caller_path[:len(caller_path) - 24])
    
    for ch in range(len(chrom_list)):
        if comp_call[ch]:
            logging.info('  signatures extracted from ' + chrom_list[ch])
        else:
            logging.info('  no signature from ' + chrom_list[ch])
        for call in comp_call[ch]:
            call[13] = [[s, read_hap[ch][s]['hap'], read_hap[ch][s]['ps'], read_hap[ch][s]['pc']] \
                        if s in read_hap[ch] else [s] for s in call[13]]
    callset = []
    for ch in range(len(chrom_list)):
        for c in range(len(comp_call[ch])):
            callset.append(dict())
            t = len(callset) - 1
            callset[t]['raw'] = comp_call[ch][c][:10]
            callset[t]['chrom'] = comp_call[ch][c][0]
            callset[t]['pos'] = int(comp_call[ch][c][1])
            callset[t]['callid'] = comp_call[ch][c][2]
            callset[t]['qual'] = comp_call[ch][c][5]
            callset[t]['filt'] = comp_call[ch][c][6]
            callset[t]['svlen'] = abs(comp_call[ch][c][10])
            callset[t]['svtype'] = comp_call[ch][c][11]
            callset[t]['svreadinfo'] = comp_call[ch][c][13]
            callset[t]['svread'] = comp_call[ch][c][12]
            callset[t]['callgt'] = comp_call[ch][c][14]
            callset[t]['refread'] = comp_call[ch][c][15]
    return callset

def get_phase_info(call, ps_sr = 2, ps_num = 1, oneps_set = ''):
    info = dict()
    info['hap1'] = info['hap2'] = info['tothap'] = info['hap1_totsc'] = info['allhap'] \
    = info['hap2_totsc'] = info['hap1_avgsc'] = info['hap2_avgsc'] = info['ps'] = info['hap0'] = 0
    if ps_num == 1:
        for read in call['svreadinfo']:
            if len(read) > 1 and read[3] <= 8100:
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
            if len(read) > 1 and read[3] <= 8100:
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

def predict_hp(call, ps_num, oneps_set):
    f = get_phase_info(call, 2, ps_num, oneps_set)
    pred = 0
    if ps_num == 0:
        if f['sv_ratio'] == 1 and f['sv_num'] >= 4:
            pred = 3
    elif ps_num == 2:
        if f['sv_ratio'] >= 0.72:
            if f['hap_avgsc_diff'] <= 1369.50:
                if f['sv_num'] >= 3:
                    pred = 3
            else:
                if f['hap0'] >= 6:
                    pred = 3
    elif ps_num == 1:
        if (f['hap1'] == 0 and f['hap2'] == 0) or f['sv_num'] >= 20:
            pred = 0
        if f['onehap_totsc'] != 0:
            if f['sv_ratio'] <= 0.24:
                pred = 0
            elif f['sv_ratio'] <= 0.9:
                if f['hapread_ratio'] <= 0.75 and f['hap_avgsc_diff'] <= 2400 or f['hapread_ratio'] > 0.75:
                    pred = 1 if f['hap1_avgsc'] > 0 else 2
            else:
                if f['hapread_ratio'] <= 0.75 and f['hap_avgsc_diff'] <= 2400 or f['hapread_ratio'] > 0.75:
                    pred = 3
        if f['onehap_totsc'] == 0:
            if f['sv_ratio'] <= 0.3:
                pred = 0
            elif f['sv_ratio'] <= 0.45:
                if f['ref_num'] > 10:
                    pred = 0
                else:
                    pred = 1 if f['hap1_totsc'] > f['hap2_totsc'] else 2
            elif f['sv_ratio'] <= 0.75:
                if f['totsc_ratio'] <= 9.72:
                    pred = 3
                else:
                    pred = 1 if f['hap1_totsc'] > f['hap2_totsc'] else 2
            else:
                pred = 3
    return pred, f['ps']

def generate_phased_callset(vcf_path, sam_home, svlen_thres, suppread_thres, thread, include_all_ctgs):
    callstat = generate_callinfo(vcf_path, read_hap_bam(sam_home, thread, include_all_ctgs), include_all_ctgs) # sam_home = home + '/snp_phasing/' -13
    logging.info('integrate read weight information')
    chrom_list = init_chrom_list(include_all_ctgs, vcf_path[:len(vcf_path) - 24])
    callset = [c for c in callstat if c['svlen'] >= svlen_thres and \
               c['svread'] >= suppread_thres and c['callgt'] not in ['./.']]
    callset_dict = dict()
    callset_dict[1] = [c for c in callset if len(set([s[2] for s in c['svreadinfo'] if len(s) > 1])) == 1]
    callset_dict[2] = [c for c in callset if len(set([s[2] for s in c['svreadinfo'] if len(s) > 1])) > 1]
    callset_dict[0] = [c for c in callset if len(set([s[2] for s in c['svreadinfo'] if len(s) > 1])) == 0]
    oneps_set = [set() for ch in range(len(chrom_list))]
    logging.info('calculate read weight statistics')
    for ch in range(len(chrom_list)):
        callset_ch = [c for c in callset_dict[1] if c['chrom'] in ['chr' + chrom_list[ch], chrom_list[ch]]]
        for c in callset_ch:
            for read in c['svreadinfo']:
                if len(read) > 1 and read[3] <= 8100:
                    oneps_set[ch].add(read[2])
                    break
    phased_callset = []
    logging.info('predict SV haplotypes in the callset')
    for ch in range(len(chrom_list)):
        for ps_num in range(3):
            callset_ch = [c for c in callset_dict[ps_num] if c['chrom'] in ['chr' + chrom_list[ch], chrom_list[ch]]]
            if not oneps_set[ch]:
                continue
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
