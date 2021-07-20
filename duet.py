import argparse
import os
import numpy as np
from operator import itemgetter
import logging
import time
import sys

def init_chrom_list():
    chrom = []
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

def parse_vcf(vcf_file, support_reads_file = ''):
    chrom_list = init_chrom_list()
    callset = read_file(vcf_file)
    call = [[] for ch in range(24)]
    for ch in range(24):
        call[ch] = [s for s in callset if s[0] in ['chr' + chrom_list[ch], chrom_list[ch]]]
        if call[ch] == []:
            continue
        info = [s[7].split(';') for s in call[ch]]
    
        svlen = [[s2 for s2 in s if 'SVLEN=' in s2] for s in info]
        svlen = ['SVLEN=0' if s == [] or s[0] == 'SVLEN=.' else s[0] for s in svlen]    
        svlen = [int(s[7:]) if '>' in s else int(s[6:]) for s in svlen]
        call[ch] = [call[ch][i] + [svlen[i]] for i in range(len(call[ch]))]

        svtype = [[s2 for s2 in s if 'SVTYPE=' in s2][0][7:] for s in info]
        call[ch] = [call[ch][i] + [svtype[i]] for i in range(len(call[ch]))]
    
        re = [[s2 for s2 in s if any(x in s2 for x in ['SUPPORT=', 'SR=', 'RE='])] for s in info]
        if re[0] != []:
            if 'SUPPORT=' in re[0][0]:
                re = [int(s[0][8:]) for s in re]
                call[ch] = [call[ch][i] + [re[i]] for i in range(len(call[ch]))]
            else:
                re = [int(s[0][3:]) for s in re]
                call[ch] = [call[ch][i] + [re[i]] for i in range(len(call[ch]))]

        rname = [[s2 for s2 in s if any(x in s2 for x in ['RNAMES=', 'READS='])] for s in info]
        if rname[0] != []:
            if 'RNAMES=' in rname[0][0]:
                rname = [s[0][7:].split(',') for s in rname]
                call[ch] = [call[ch][i] + [rname[i]] for i in range(len(call[ch]))]
            else:
                rname = [s[0][6:].split(',') for s in rname]
                call[ch] = [call[ch][i] + [rname[i]] for i in range(len(call[ch]))]
    
        if support_reads_file != '':
            raw_rname = read_file(support_reads_file)
            raw_rname = raw_rname[1:]
            nanovar_sv2read = dict()
            for s in raw_rname:
                nanovar_sv2read[s[0]] = {}
                tmp_reads = [s2[:s2.find('~')] for s2 in s[1].split(',') if len(s2) > 5]
                tmp_reads = list(dict.fromkeys(tmp_reads))
                nanovar_sv2read[s[0]] = tmp_reads
            call[ch] = [call[ch][i] + [nanovar_sv2read[call[ch][i][2]]] for i in range(len(call[ch]))]
            
        gtinfo = [s[9].split(':') for s in call[ch]]
        if len(gtinfo[0]) > 4: # gt ref var
            call[ch] = [call[ch][i] + [gtinfo[i][0]] for i in range(len(call[ch]))]
            call[ch] = [call[ch][i] + [int(gtinfo[i][1]) if gtinfo[i][1] != '.' else 0] \
                        for i in range(len(call[ch]))]
            call[ch] = [call[ch][i] + [int(gtinfo[i][2]) if gtinfo[i][2] != '.' else 0] \
                        for i in range(len(call[ch]))]
        elif len(gtinfo[0]) >= 3:
            call[ch] = [call[ch][i] + [gtinfo[i][0]] for i in range(len(call[ch]))]
            if gtinfo[0][-1].find(',') == -1:
                call[ch] = [call[ch][i] + [int(gtinfo[i][1]) if gtinfo[i][1] != '.' else 0] \
                            for i in range(len(call[ch]))]
                call[ch] = [call[ch][i] + [int(gtinfo[i][2]) if gtinfo[i][2] != '.' else 0] \
                            for i in range(len(call[ch]))]
            else:
                call[ch] = [call[ch][i] + [int(gtinfo[i][-1][:gtinfo[i][-1].find(',')]) \
                            if gtinfo[i][-1][:gtinfo[i][-1].find(',')] != '.' else 0] \
                            for i in range(len(call[ch]))]
                call[ch] = [call[ch][i] + [int(gtinfo[i][-1][gtinfo[i][-1].find(',') + 1:]) \
                            if gtinfo[i][-1][gtinfo[i][-1].find(',') + 1:] != '.' else 0] \
                            for i in range(len(call[ch]))]
    return(call)

def read_hap_bam(path):
    logging.info('extract read weight information from haplotagged alignments')
    read_hap = [[] for ch in range(24)]
    chrom_list = init_chrom_list()
    for ch in range(24):
        read_hap[ch] = dict()
        if os.path.exists(path + 'chr' + chrom_list[ch] + '.sam'):
            hap_bam_path = path + 'chr' + chrom_list[ch] + '.sam'
        else:
            hap_bam_path = path + chrom_list[ch] + '.sam'
        aln_info = [[s[0], int(s[1]), int(s[3]), int(s[4]), int(s[-3][5:]), int(s[-2][5:]), int(s[-1][5:])] \
                        for s in read_file(hap_bam_path) if 'PC:i:' in s[-2]]
        for item in aln_info:
            read_hap[ch][item[0]] = {'hap': item[-3], 'ps': item[-1], 'pc': item[-2]}
    return read_hap

def generate_callinfo(caller_path, read_hap):
    logging.info('extract read weight information from structural variation signatures')
    comp_call = parse_vcf(caller_path)
    for ch in range(24):
        for call in comp_call[ch]:
            call[13] = [[s, read_hap[ch][s]['hap'], read_hap[ch][s]['ps'], read_hap[ch][s]['pc']] \
                        if s in read_hap[ch] else [s] for s in call[13]]
    callset = []
    for ch in range(24):
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

def generate_phased_callset(vcf_path, sam_home, svlen_thres, suppread_thres):
    callstat = generate_callinfo(vcf_path, read_hap_bam(sam_home))
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
                if len(read) > 1 and read[3] <= 8100:
                    oneps_set[ch].add(read[2])
                    break
    phased_callset = []
    logging.info('predict SV haplotypes in the callset')
    for ch in range(24):
        for ps_num in range(3):
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

def print_sv(phased_callset, output_path):
    logging.info('write phased callset into .vcf file')
    f = open(output_path, 'a')
    idx = 0
    for c in phased_callset:
        idx += 1
        f.write(c['chrom'] + '\t' + str(c['pos']) + '\tDuet.' + str(idx) + '\tN\t<' + c['alt'] + \
                '>\t.\tPASS\tSVLEN=' + str(c['svlen']) + '\tHP:PS\t' + c['hp'] + ':' + str(c['ps']) + '\n')
    f.close()
    
def print_sv_header(vcf_path, output_path):
    alt_str = """##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">
##ALT=<ID=DEL,Description="Deletion relative to the reference">
"""
    filter_str = '##FILTER=<ID=PASS,Description="SV calls passed phasing criterion">\n'
    info_str = '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Estimated length of the variant">\n'
    format_str = """##FORMAT=<ID=HP,Number=1,Type=String,Description="Haplotype of the SV call">
##FORMAT=<ID=PS,Number=1,Type=String,Description="Phase set which the SV call belongs to">"""
    fileformat_str = '##fileformat=VCFv4.2\n'
    source_str = '##source=Duet\n'
    header_str = fileformat_str + source_str + alt_str + filter_str + info_str + format_str + '\n'
    vcf_str = read_file(vcf_path)
    chrom_list = init_chrom_list()
    contig_str = ['' for ch in range(24)]
    for ch in range(24):
        for l in vcf_str:
            if '##contig=<ID=chr' + chrom_list[ch] + ',' in l[0] or '##contig=<ID=' + chrom_list[ch] + ',' in l[0]:
                header_str += l[0] + '\n'
                continue
    header_str += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVALUE\n'
    f = open(output_path, 'w')
    f.write(header_str)
    f.close()
    
def snp_phasing(home, ref_path, aln_path, thread):
    lines = '*************************'
    logging.info(lines + ' SNP PHASING STARTED ' + lines)
    starttime = time.time()
    snp_calling_home = home + '/snp_calling/'
    snp_phasing_home = home + '/snp_phasing/'
    os.system('mkdir ' + snp_phasing_home)
    pileup_vcf_path = snp_calling_home + 'pileup.vcf'
    os.system('gunzip -c ' + pileup_vcf_path + '.gz > ' + pileup_vcf_path)
    chrom_list = init_chrom_list()
    vcf_file = read_file(pileup_vcf_path)
    vcf_header = []
    vcf_call = []
    for l in vcf_file:
        if l[0][0] == '#':
            vcf_header.append(l)
        else:
            vcf_call.append(l)
    for ch in range(24):
        chrom_vcf = [c for c in vcf_call if c[0] in ['chr' + chrom_list[ch], chrom_list[ch]] and c[6] != 'RefCall']
        fpath = snp_phasing_home + chrom_vcf[0][0] + '.vcf'
        f = open(fpath, 'w')
        for l in vcf_header:
            if l[0] == '#CHROM':
                st = '\t'.join(l) + '\n'
            else:
                st = ' '.join(l) + '\n'
            f.write(st)
        for l in chrom_vcf:
            st = '\t'.join(l) + '\n'
            f.write(st)
        f.close()
        # os.system('bcftools sort ' + fpath)
    f = open(snp_phasing_home + 'parallel_wh.sh', 'w')
    if chrom_vcf[0][0] in chrom_list:
        f.write('CHR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)\n')
    else:
        f.write('CHR=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)\n')
    f.write('parallel -j' + str(thread) + ' "whatshap phase -o ' + snp_phasing_home + 'phased_{1}.vcf.gz -r ' + \
            ref_path + ' --chromosome {1} --distrust-genotypes --ignore-read-groups ' + \
            snp_phasing_home + '{1}.vcf ' + aln_path + '" ::: ${CHR[@]}\n')
    f.write('parallel -j' + str(thread) + ' "tabix -f -p vcf ' + snp_phasing_home + \
            'phased_{1}.vcf.gz" ::: ${CHR[@]}\n')
    f.write('parallel -j' + str(thread) + ' "whatshap haplotag -o ' + snp_phasing_home + '{1}.bam -r ' + \
            ref_path + ' --regions {1} --ignore-read-groups --tag-supplementary ' + \
            snp_phasing_home + 'phased_{1}.vcf.gz ' + aln_path + '" ::: ${CHR[@]}\n')
    f.write('parallel -j' + str(thread) + ' "samtools index -@' + str(thread) + ' ' + \
            snp_phasing_home + '{1}.bam" ::: ${CHR[@]}\n')
    f.write('parallel -j' + str(thread) + ' "samtools view -h -@' + str(thread) + ' ' + snp_phasing_home + \
            '{1}.bam > ' + snp_phasing_home + '{1}.sam" ::: ${CHR[@]}\n')
    f.close()
    os.system('chmod a+x ' + snp_phasing_home + 'parallel_wh.sh')
    os.system(snp_phasing_home + 'parallel_wh.sh')
    logging.info(lines + ' SNP PHASING COMPLETED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)
    
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
    
def sv_phasing(home, svlen_thres, suppread_thres):
    lines = '*************************'
    logging.info(lines + ' SV PHASING STARTED ' + lines)
    starttime = time.time()
    sv_calling_path = home + '/sv_calling/variants.vcf'
    sv_phasing_path = home + '/phased_sv.vcf'
    snp_phasing_home = home + '/snp_phasing/'
    logging.info('create output .vcf file')
    print_sv_header(sv_calling_path, sv_phasing_path)
    print_sv(generate_phased_callset(sv_calling_path, snp_phasing_home, svlen_thres, suppread_thres), \
             sv_phasing_path)
    logging.info(lines + ' SV PHASING COMPLETED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)

def snp_calling(home, ref_path, aln_path, maf):
    lines = '*************************'
    logging.info(lines + ' SNP CALLING STARTED ' + lines)
    starttime = time.time()
    snp_calling_home = home + '/snp_calling/'
    os.system('mkdir ' + snp_calling_home)
    os.system('run_clair3.sh -b ' + aln_path + ' -f ' + ref_path + \
              ' -m "${CONDA_PREFIX}/bin/models/ont" -t 40 -p ont -o ' + snp_calling_home + \
              ' --snp_min_af=' + str(maf) + ' --pileup_only --call_snp_only')
    logging.info(lines + ' SNP CALLING COMPLETED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)    

def set_logging(home):
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt = '%H:%M:%S')
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)
    fileHandler = logging.FileHandler(home + '/run_duet.log', mode="w")
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description = 'Duet: a fast and accurate germline variant caller for structural variation calling and haplotyping from low-depth nanopore sequencing reads.')
    parser.add_argument('-t', '--thread', type = int, default = 40, 
        help = 'number of threads to use [%(default)s]')
    parser.add_argument('-m', '--min_allele_frequency', type = float, default = 0.25, 
        help = 'minimum allele frequency required to call a candidate SNP [%(default)s]')
    parser.add_argument('-c', '--cluster_max_distance', type = float, default = 0.9, 
        help = 'maximum span-position distance between SV marks in a cluster to call a SV candidates [%(default)s]')
    parser.add_argument('-s', '--sv_min_size', type = int, default = 50, 
        help = 'minimum SV size to be reported [%(default)s]')
    parser.add_argument('-r', '--min_support_read', type = int, default = 2, 
        help = 'minimum number of reads that support a SV to be reported [%(default)s]')
    parser.add_argument('BAM', type = str, 
	    help = 'sorted alignment file in .bam format (along with .bai file in the same directory)')
    parser.add_argument('REFERENCE', type = str, 
	    help = 'indexed reference genome in .fasta format (along with .fai file in the same directory)')
    parser.add_argument('OUTPUT', type = str, 
        help = 'working and output directory (existing files in the directory will be overwritten)')
    args = parser.parse_args()
    return args

def run(argv):
    args = parse_args(argv)
    if not os.path.exists(args.OUTPUT):
        os.makedirs(args.OUTPUT)
    set_logging(args.OUTPUT)
    lines = '*************************'
    starttime = time.time()
    logging.info(lines + ' DUET STARTED ' + lines)
    snp_calling(args.OUTPUT, args.REFERENCE, args.BAM, args.min_allele_frequency)
    sv_calling(args.OUTPUT, args.REFERENCE, args.BAM, args.cluster_max_distance, args.sv_min_size)
    snp_phasing(args.OUTPUT, args.REFERENCE, args.BAM, args.thread)
    sv_phasing(args.OUTPUT, args.sv_min_size, args.min_support_read)
    logging.info(lines + ' DUET FINISHED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)
    logging.info('OUTPUT .VCF FILE AT ' + args.OUTPUT + 'phased_sv.vcf')

if __name__ == '__main__':
    run(sys.argv[1:])