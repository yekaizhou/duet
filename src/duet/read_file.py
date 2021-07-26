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