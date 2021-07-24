from read_file import init_chrom_list, read_file
import logging

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