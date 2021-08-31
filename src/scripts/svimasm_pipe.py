import os
import shlex
import subprocess
import sys
import argparse
import shutil

def init_chrom_list():
    chrom_list = []
    for i in range(1, 23):
        chrom_list.append(str(i))
    chrom_list.append('X')
    chrom_list.append('Y')
    return(chrom_list)

def parse_args(argv):
    parser = argparse.ArgumentParser(description = 'SV phasing pipeline built on svim-asm')
    parser.add_argument('working_dir', type = str,
            help = 'working and output directory')
    parser.add_argument('hap_bam', type = str,
            help = 'haplotagged bam file')
    parser.add_argument('ref', type = str,
            help = 'indexed reference genome in .fasta format (along with .fai file in the same directory)')
    args = parser.parse_args()
    return args

def main(argv):
    args = parse_args(argv)

    thread = 40
    home = args.working_dir
    ref_home = args.ref
    hap_bam_path = args.hap_bam

    tmp_sam_hp1_path = home + '/tmp_hp1.sam'
    tmp_sam_hp2_path = home + '/tmp_hp2.sam'
    tmp_dipcall_path = home + '/variants.vcf'
    dipcallset_path = home + '/phased_sv.vcf'

    aln_header = subprocess.check_output(shlex.split('samtools view -H -@' + str(thread) + ' ' + hap_bam_path)).decode('ascii').split('\n')[:-1]
    hap_bam_raw = subprocess.check_output(shlex.split('samtools view -@' + str(thread) + ' ' + hap_bam_path)).decode('ascii').split('\n')[:-1]
    hap_bam_raw = [s.split() for s in hap_bam_raw]
    hap_bam_raw = [s for s in hap_bam_raw if 'PC:i:' in s[-2]]
    chrom_list = init_chrom_list()
    hap_bam = [[s for s in hap_bam_raw if s[2] in [chrom_list[ch], 'chr' + chrom_list[ch]]] for ch in range(24)]

    tot_vcf_str = ''

    for ch in range(24):
        ps_set = set([s[-1] for s in hap_bam[ch]])
        for ps in ps_set:
            bam_ps = [s for s in hap_bam[ch] if s[-1] == ps]
            bam_ps_hp1 = [s for s in bam_ps if s[-3] == 'HP:i:1']
            bam_ps_hp2 = [s for s in bam_ps if s[-3] == 'HP:i:2']
            sam_hp1 = open(tmp_sam_hp1_path, 'w')
            sam_hp2 = open(tmp_sam_hp2_path, 'w')

            for s in aln_header:
                sam_hp1.write(s + '\n')
                sam_hp2.write(s + '\n')
    
            for s in bam_ps_hp1:
                aln_str = '\t'.join(s) + '\n'
                sam_hp1.write(aln_str)
    
            for s in bam_ps_hp2:
                aln_str = '\t'.join(s) + '\n'
                sam_hp2.write(aln_str)
    
            sam_hp1.close()
            sam_hp2.close()
            
            os.system('samtools fasta -@' + str(thread) + ' ' + tmp_sam_hp1_path + ' > ' + tmp_sam_hp1_path[:-4] + '.fa')
            os.system('samtools fasta -@' + str(thread) + ' ' + tmp_sam_hp2_path + ' > ' + tmp_sam_hp2_path[:-4] + '.fa')
            os.system('flye --nano-raw ' + tmp_sam_hp1_path[:-4] + '.fa --out-dir ' + home + '/flye_hp1 -t ' + str(thread))
            os.system('flye --nano-raw ' + tmp_sam_hp2_path[:-4] + '.fa --out-dir ' + home + '/flye_hp2 -t ' + str(thread))
            os.system('minimap2 -a -x asm5 --cs -r2k -t ' + str(thread) + ' ' + ref_home + ' ' + home + '/flye_hp1/assembly.fasta > ' + home + '/flye_hp1/assembly.sam')
            os.system('minimap2 -a -x asm5 --cs -r2k -t ' + str(thread) + ' ' + ref_home + ' ' + home + '/flye_hp2/assembly.fasta > ' + home + '/flye_hp2/assembly.sam')
            os.system('samtools sort -m4G -@ ' + str(thread) + ' -o ' + home + '/flye_hp1/assembly.bam ' + home + '/flye_hp1/assembly.sam')
            os.system('samtools sort -m4G -@ ' + str(thread) + ' -o ' + home + '/flye_hp2/assembly.bam ' + home + '/flye_hp2/assembly.sam')
            os.system('samtools index -@ ' + str(thread) + ' ' + home + '/flye_hp1/assembly.bam')
            os.system('samtools index -@ ' + str(thread) + ' ' + home + '/flye_hp2/assembly.bam')
            os.system('svim-asm diploid --min_sv_size 50 ' + home + ' ' + home + '/flye_hp1/assembly.bam ' + home + '/flye_hp2/assembly.bam ' + ref_home)
            os.system('rm -rf ' + home + '/flye_hp1 ' + home + '/flye_hp2')
            if not os.path.exists(tmp_dipcall_path):
                continue 
            vcf_header = subprocess.check_output(shlex.split('bcftools view -h --threads ' + str(thread) + ' ' + tmp_dipcall_path)).decode('ascii').split('\n')[:-1]
            tmp_vcf_str = subprocess.check_output(shlex.split('bcftools view -H --threads ' + str(thread) + ' ' + tmp_dipcall_path)).decode('ascii').split('\n')[:-1]
            for s in tmp_vcf_str:
                tot_vcf_str += s[:-2] + '|' + s[-1] + ':' + str(ps[5:]) + '\n'
            
    dip_vcf = open(dipcallset_path, 'w')
    for s in vcf_header:
        dip_vcf.write(s + '\n')
    dip_vcf.write(tot_vcf_str)
    dip_vcf.close()
    
if __name__ == '__main__':
    main(sys.argv[1:])
