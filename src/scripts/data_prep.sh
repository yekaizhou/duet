# make sure there are minimap2, samtools, and bcftools etc in current environment

# enter a folder for data download, storage and process
mkdir duet_bench && cd duet_bench

# download and process hg38 reference genome used for HGSVC2 SV characterization
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz
gunzip -d -k hg38.no_alt.fa.gz
sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' hg38.no_alt.fa
samtools faidx hg38.no_alt.fa

# download and process HG001, HG002, and HG00733 phased SV truth sets generated by HGSVC2
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_sv_insdel_sym.vcf.gz
gunzip -d -k variants_freeze4_sv_insdel_sym.vcf.gz
sed -i 's/\./PASS/2' variants_freeze4_sv_insdel_sym.vcf
grep -v '_GL\|_KI\|chrM' variants_freeze4_sv_insdel_sym.vcf > variants_freeze4_sv_insdel_sym_mainChr.vcf
bcftools view -p -c1 -Ov -s NA12878 -o NA12878_HGSVC2.vcf variants_freeze4_sv_insdel_sym_mainChr.vcf
bcftools view -p -c1 -Ov -s NA24385 -o NA24385_HGSVC2.vcf variants_freeze4_sv_insdel_sym_mainChr.vcf
bcftools view -p -c1 -Ov -s HG00733 -o HG00733_HGSVC2.vcf variants_freeze4_sv_insdel_sym_mainChr.vcf

# download and process HG002 CMRG medically-relevant SV truth set
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.bed
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz
gunzip HG002_GRCh38_CMRG_SV_v1.00.vcf.gz
mv HG002_GRCh38_CMRG_SV_v1.00.vcf HG002_CMRG.vcf
mv HG002_GRCh38_CMRG_SV_v1.00.bed HG002_CMRG.bed

# download HG001, HG002, and HG00733 ONT reads called by Guppy4.2.2
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG001/nanopore/Guppy_4.2.2/HG001_NBT2018_Guppy_4.2.2.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/Guppy_4.2.2/HG002_GIAB_PromethION_Guppy_4.2.2_prom.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG00733/nanopore/Guppy_4.2.2/HG00733_1_Guppy_4.2.2_prom.fastq.gz

minimap2 -ax map-ont hg38.no_alt.fa HG001_NBT2018_Guppy_4.2.2.fastq.gz -t 40 -o HG001_NBT2018_Guppy_4.2.2.sam
minimap2 -ax map-ont hg38.no_alt.fa HG002_GIAB_PromethION_Guppy_4.2.2_prom.fastq.gz -t 40 -o HG002_GIAB_PromethION_Guppy_4.2.2_prom.sam
minimap2 -ax map-ont hg38.no_alt.fa HG00733_1_Guppy_4.2.2_prom.fastq.gz -t 40 -o HG00733_1_Guppy_4.2.2_prom.sam

samtools view -buS -@40 HG001_NBT2018_Guppy_4.2.2.sam | samtools sort -@40 -O bam -T ./ - > HG001_NBT2018_Guppy_4.2.2.bam
samtools view -buS -@40 HG002_GIAB_PromethION_Guppy_4.2.2_prom.sam | samtools sort -@40 -O bam -T ./ - > HG002_GIAB_PromethION_Guppy_4.2.2_prom.bam
samtools view -buS -@40 HG00733_1_Guppy_4.2.2_prom.sam | samtools sort -@40 -O bam -T ./ - > HG00733_1_Guppy_4.2.2_prom.bam

samtools view -bS -s 0.23992754 -@40 HG001_NBT2018_Guppy_4.2.2.bam > NA12878_8X.bam
samtools view -bS -s 0.52230572 -@40 HG002_GIAB_PromethION_Guppy_4.2.2_prom.bam > NA24385_8X.bam
samtools view -bS -s 0.42983715 -@40 HG00733_1_Guppy_4.2.2_prom.bam > HG00733_8X.bam

samtools index -@40 NA12878_8X.bam
samtools index -@40 NA24385_8X.bam
samtools index -@40 HG00733_8X.bam

# full depth version
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG001/nanopore/Guppy_4.2.2/HG001_Circulomics_Guppy_4.2.2.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/Guppy_4.2.2/HG002_GIAB_MinION_GridION_Guppy_4.2.2.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/Guppy_4.2.2/HG002_GIAB_PromethION_Guppy_4.2.2_prom.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG00733/nanopore/Guppy_4.2.2/HG00733_1_Guppy_4.2.2_prom.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG00733/nanopore/Guppy_4.2.2/HG00733_2_Guppy_4.2.2_prom.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG00733/nanopore/Guppy_4.2.2/HG00733_3_Guppy_4.2.2_prom.fastq.gz

minimap2 -ax map-ont hg38.no_alt.fa HG001_Circulomics_Guppy_4.2.2.fastq.gz -t 40 -o HG001.sam
minimap2 -ax map-ont hg38.no_alt.fa HG002_GIAB_* -t 60 -o HG002.sam
minimap2 -ax map-ont hg38.no_alt.fa HG00733_1_Guppy_4.2.2_prom.fastq.gz HG00733_2_Guppy_4.2.2_prom.fastq.gz -t 40 -o HG00733.sam





# results in files for further utilization:
#   Reference: hg38.no_alt.fa
#   Truth SV callset: NA12878_HGSVC2.vcf NA24385_HGSVC2.vcf HG00733_HGSVC2.vcf HG002_CMRG.vcf HG002_CMRG.bed
#   Alignment: NA12878_8X.bam NA24385_8X.bam HG00733_8X.bam
