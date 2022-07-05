from os import read
import sys
import argparse
import numpy as np
import shlex
import subprocess
import numpy as np
from operator import itemgetter

def init_chrom_list():
    chrom_list = []
    for i in range(1, 23):
        chrom_list.append(str(i))
    chrom_list.append('X')
    chrom_list.append('Y')
    return(chrom_list)

def read_file(vcf_path):
    with open(vcf_path, 'r') as file:
        raw = file.readlines()
    raw = [s.strip() for s in raw]
    raw = [s.split() for s in raw]
    return raw

def parse_bed(path):
    bed_raw = read_file(path)
    chrom_list = init_chrom_list()
    comp_range = [[] for ch in range(24)]
    for ch in range(24):
        bed_raw_ch = [s for s in bed_raw if s[0] == 'chr' + chrom_list[ch]]
        for s in bed_raw_ch:
            comp_range[ch].append(range(int(s[1]), int(s[2]) + 1))
    return comp_range

def parse_vcf(vcf_path, skip_phasing, bed_path = ''):
    if bed_path != '':
        comp_range = parse_bed(bed_path)
    chrom_list = init_chrom_list()
    sv_set = read_file(vcf_path)
    sv_set = [s for s in sv_set if s[0][0] != '#']
    info = []
    for s in sv_set:
        if s[0][3:] not in chrom_list:
            continue
        if 'SVLEN=.' in s[7]:
            continue
        if 'INS' not in s[7] and 'DEL' not in s[7] and 'DUP' not in s[7] and 'INS' not in s[4] and 'DEL' not in s[4] and 'DUP' not in s[4]:
            continue
        info.append(dict())
        t = len(info) - 1
        info[t]['chr'] = s[0]
        info[t]['pos'] = int(s[1])
        info[t]['id'] = s[2] + s[0] + s[1]
        info[t]['hp'] = s[-1][:3]
        #if t == 0:
        #    info[t]['ps'] = s[0]
        if info[t]['hp'][0] == '.':
            info[t]['hp'] = '0' + info[t]['hp'][1:]
        if info[t]['hp'][2] == '.':
            info[t]['hp'] = info[t]['hp'][:2] + '0'
        if info[t]['hp'][1] == '/':
            if not skip_phasing:
                if info[t]['hp'] == '1/1':
                    info[t]['ps'] = s[0] + '_' + s[-1][s[-1].rfind(':'):]
                    #info[t]['ps'] = info[t - 1]['ps']
                else:
                    info.pop()
                    continue
            info[t]['hp'] = info[t]['hp'][0] + '|' + info[t]['hp'][2]
            info[t]['ps'] = s[0]
        else:
            info[t]['ps'] = s[0] + '_' + s[-1][s[-1].rfind(':'):] # baseinfo ps will be only chrN_
        sv_info = s[7].split(';')
        if 'SVLEN' in s[7]:
            info[t]['len'] = abs(int([s for s in sv_info if 'SVLEN' in s][0][6:])) if 'SVLEN=>' not in s[7] else abs(int([s for s in sv_info if 'SVLEN' in s][0][7:]))
            info[t]['type'] = s[4][1:-1] if s[4] in ['<INS>', '<DEL>', '<DUP:TANDEM>', '<DUP:INT>'] else [s for s in sv_info if 'SVTYPE' in s][0][7:]
            if 'DUP' in info[t]['type']:
                info[t]['type'] = 'INS'
        else:
            l = len(s[3]) - len(s[4])
            if l > 0:
                info[t]['len'] = l
                info[t]['type'] = 'DEL'
            if l < 0:
                info[t]['len'] = -l
                info[t]['type'] = 'INS'
        if bed_path == '':
            flag = True
        else:
            flag = False
            for bed_rg in comp_range[chrom_list.index(s[0][3:])]:
                if int(s[1]) in bed_rg:
                    flag = True
                    break
        if not flag or info[t]['len'] < 50 or info[t]['hp'] == '0|0':
            info.pop()
    return info

def evaluation(baseinfo, callinfo, threshold_tp_range, ratio):
    chrom_list = init_chrom_list()
    call_tp, call_tp_gt, call_tp_hp, base_tp, base_tp_gt, base_tp_hp = set(), set(), set(), set(), set(), set()
    avg_sv_num = len(callinfo) / len(set([s['ps'] for s in callinfo]))
    for ch in range(24):
        base_ch_type = dict()
        for svtype in ['INS', 'DEL']:
            base_ch_type[svtype] = sorted([s for s in baseinfo if s['chr'] == 'chr' + chrom_list[ch] and s['type'] == svtype], key = itemgetter('pos'))
        call_ch = [s for s in callinfo if s['chr'] == 'chr' + chrom_list[ch]]
        ps_set = set([s['ps'] for s in call_ch])
        for ps in ps_set:
            call_ch_ps = [s for s in call_ch if s['ps'] == ps]
            tmp1_call_tp_hp, tmp1_base_tp_hp, tmp2_call_tp_hp, tmp2_base_tp_hp = set(), set(), set(), set()
            for svtype in ['INS', 'DEL']:
                call = [s for s in call_ch_ps if s['type'] == svtype]
                base = base_ch_type[svtype]
                if not call:
                    continue
                idx_list = np.searchsorted([s['pos'] for s in base], [s['pos'] for s in call])
                for call_idx in range(len(idx_list)):
                    if idx_list[call_idx] == len(base):
                        base_idx = idx_list[call_idx] - 1
                    elif idx_list[call_idx] > 0 and abs(call[call_idx]['pos'] - base[idx_list[call_idx]]['pos']) > \
                        abs(call[call_idx]['pos'] - base[idx_list[call_idx] - 1]['pos']):
                        base_idx = idx_list[call_idx] - 1
                    else:
                        base_idx = idx_list[call_idx]
                    if abs(call[call_idx]['pos'] - base[base_idx]['pos']) <= threshold_tp_range and \
                       (min(call[call_idx]['len'], base[base_idx]['len']) / max(base[base_idx]['len'], call[call_idx]['len']) >= ratio):
                        call_tp.add(call[call_idx]['id'])
                        base_tp.add(base[base_idx]['id'])
                        if call[call_idx]['hp'] in ['1|0', '0|1'] and base[base_idx]['hp'] in ['1|0', '0|1'] or \
                           call[call_idx]['hp'] == base[base_idx]['hp'] == '1|1':
                            call_tp_gt.add(call[call_idx]['id'])
                            base_tp_gt.add(base[base_idx]['id'])
                        if call[call_idx]['hp'] == base[base_idx]['hp']:
                            tmp1_call_tp_hp.add(call[call_idx]['id'])
                            tmp1_base_tp_hp.add(base[base_idx]['id'])
                        if call[call_idx]['hp'] == base[base_idx]['hp'] == '1|1' or \
                           call[call_idx]['hp'] == '0|1' and base[base_idx]['hp'] == '1|0' or \
                           call[call_idx]['hp'] == '1|0' and base[base_idx]['hp'] == '0|1':
                            tmp2_call_tp_hp.add(call[call_idx]['id'])
                            tmp2_base_tp_hp.add(base[base_idx]['id'])
                            
            if len(tmp1_call_tp_hp) + len(tmp1_base_tp_hp) > len(tmp2_call_tp_hp) + len(tmp2_base_tp_hp):
                call_tp_hp = call_tp_hp.union(tmp1_call_tp_hp)
                base_tp_hp = base_tp_hp.union(tmp1_base_tp_hp)
            else:
                call_tp_hp = call_tp_hp.union(tmp2_call_tp_hp)
                base_tp_hp = base_tp_hp.union(tmp2_base_tp_hp)
            
    p = len(call_tp) / len(callinfo)
    r = len(base_tp) / len(baseinfo)
    f1 = 2 * p * r / (p + r)
    p_gt = len(call_tp_gt) / len(callinfo)
    r_gt = len(base_tp_gt) / len(baseinfo)
    f1_gt = 2 * p_gt * r_gt / (p_gt + r_gt)
    p_hp = len(call_tp_hp) / len(callinfo)
    r_hp = len(base_tp_hp) / len(baseinfo)
    f1_hp = 2 * p_hp * r_hp / (p_hp + r_hp)
    return avg_sv_num, p, r, f1, p_gt, r_gt, f1_gt, p_hp, r_hp, f1_hp

def parse_args(argv):
    parser = argparse.ArgumentParser(description = 'evaluate SV calling, genotyping and phasing performance')
    parser.add_argument('callset', type = str,
            help = 'phased SV callset in .vcf format')
    parser.add_argument('truthset', type = str,
            help = 'phased SV truthset in .vcf format')
    parser.add_argument('-r', '--refdist', type = int, default = 1000,
            help = 'maximum distance comparison calls must be within from base call')
    parser.add_argument('-p', '--pctsim', type = float, default = 0,
            help = 'edit distance ratio between the REF/ALT haplotype sequences of base and comparison call')
    parser.add_argument('-b', '--bed_file', type = str, 
            help = 'optional .bed file to confine benchmark regions')
    parser.add_argument('--skip_phasing', action = 'store_true',
            help = 'only benchmark on SV calling and genotyping [%(default)s]')
    args = parser.parse_args()
    return args

def main(argv):
    args = parse_args(argv)
    if not args.bed_file:
        avg_sv_num, p, r, f1, p_gt, r_gt, f1_gt, p_hp, r_hp, f1_hp = evaluation(parse_vcf(args.truthset, args.skip_phasing, ''), parse_vcf(args.callset, args.skip_phasing, ''), args.refdist, args.pctsim)
    else:
        avg_sv_num, p, r, f1, p_gt, r_gt, f1_gt, p_hp, r_hp, f1_hp = evaluation(parse_vcf(args.truthset, args.skip_phasing, args.bed_file), parse_vcf(args.callset, args.skip_phasing, args.bed_file), args.refdist, args.pctsim)
    if not args.skip_phasing:
        print('Average SV number per phase set is', avg_sv_num)
    print('The precision, recall and F1 score of SV calling are', p, r, f1)
    print('The precision, recall and F1 score of SV genotyping are', p_gt, r_gt, f1_gt)
    if not args.skip_phasing:
        print('The precision, recall and F1 score of SV phasing are', p_hp, r_hp, f1_hp)
    
if __name__ == '__main__':
    main(sys.argv[1:])
