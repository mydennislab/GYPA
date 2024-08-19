
# GYPA/GYPB SVs in chimpanzee

## 1. Hg38 as reference

We genotyped two candidate SNPs and asses if they were in LD with a SV. The two candidate SNPs we identify in GYPA are chr4:145039806 and chr4:145040845 in Hg19.

In GRCh38, the coordinates are:
- chr4:145039806 -> chr4:144118653 
- chr4:145040845 -> chr4:144119692

### 1.1 Assemblies

```bash
cd /share/dennislab/users/dcsoto/other/GYPA/assemblies
conda activate assembly
```

Downloading assemblies:
```bash
wget https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/mPanTro3.pri.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/mPanTro3.alt.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/mPanPan1.pri.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/mPanPan1.alt.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/mGorGor1.pri.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/mGorGor1.alt.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/mPonAbe1.pri.cur.20231205.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/mPonAbe1.alt.cur.20231205.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/mPonPyg2.pri.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/mPonPyg2.alt.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/mSymSyn1.pri.cur.20231205.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/mSymSyn1.alt.cur.20231205.fasta.gz

for file in *.20231122.fasta.gz; do gunzip $file; done
for file in *.20231205.fasta.gz; do gunzip $file; done

for file in *.20231122.fasta; do samtools faidx $file; done
for file in *.20231205.fasta; do samtools faidx $file; done
```

#### 1.1.1 Alignments

First, we compared chimpanzee assemblies to Hg38 using alignments.

```bash
# AG18354 (ToLID: mPanTro3) - stringent mapping
wget https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/assembly_verkko_1.1-hic-freeze/mPanTro3.hap1.20221111.fasta.gz
wget https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/assembly_verkko_1.1-hic-freeze/mPanTro3.hap1.20221111.fasta.gz.fai
wget https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/assembly_verkko_1.1-hic-freeze/mPanTro3.hap2.20221111.fasta.gz
wget https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/assembly_verkko_1.1-hic-freeze/mPanTro3.hap2.20221111.fasta.gz.fai

minimap2 -t 64 -a --eqx --cs -x asm5 --secondary=no -s 25000 -K 8G ../references/hg38.noalt.fa mPanTro3.hap1.20221111.fasta.gz | samtools view -F 4 -b - | samtools sort > AG18354.hap1.hg38.bam
minimap2 -t 64 -a --eqx --cs -x asm5 --secondary=no -s 25000 -K 8G ../references/hg38.noalt.fa mPanTro3.hap2.20221111.fasta.gz | samtools view -F 4 -b - | samtools sort > AG18354.hap2.hg38.bam

samtools index AG18354.hap1.hg38.bam
samtools index AG18354.hap2.hg38.bam

samtools view -Sb AG18354.hap1.hg38.bam chr4:143723743-144212365 > AG18354.hap1.hg38.GYPA_GYPB.bam
samtools view -Sb AG18354.hap2.hg38.bam chr4:143723743-144212365 > AG18354.hap2.hg38.GYPA_GYPB.bam 

samtools index AG18354.hap1.hg38.GYPA_GYPB.bam 
samtools index AG18354.hap2.hg38.GYPA_GYPB.bam

# AG18354 (ToLID: mPanTro3) - permissive mapping
minimap2 -ax asm5 -t 64 --eqx ../references/hg38.noalt.fa mPanTro3.hap1.20221111.fasta.gz | samtools sort -O BAM - > AG18354.hap1.hg38.plotsr.bam
minimap2 -ax asm5 -t 64 --eqx ../references/hg38.noalt.fa mPanTro3.hap2.20221111.fasta.gz | samtools sort -O BAM - > AG18354.hap2.hg38.plotsr.bam

samtools index AG18354.hap1.hg38.plotsr.bam
samtools index AG18354.hap2.hg38.plotsr.bam

samtools view -Sb AG18354.hap1.hg38.plotsr.bam chr4:143723743-144212365 > AG18354.hap1.hg38.plotsr.GYPA_GYPB.bam
samtools view -Sb AG18354.hap2.hg38.plotsr.bam chr4:143723743-144212365 > AG18354.hap2.hg38.plotsr.GYPA_GYPB.bam 

samtools index AG18354.hap1.hg38.plotsr.GYPA_GYPB.bam
samtools index AG18354.hap2.hg38.plotsr.GYPA_GYPB.bam

# AG18359
minimap2 -t 64 -a --eqx --cs -x asm5 --secondary=no -s 25000 -K 8G ../references/hg38.noalt.fa /share/dennislab/users/dcsoto/ms_asm/1_asm_shasta/ShastaRun_Jan2023/Assembly.fasta | samtools view -F 4 -b - | samtools sort > AG18359.asm.hg38.bam
samtools index AG18359.asm.hg38.bam
samtools view -Sb AG18359.asm.hg38.bam chr4:143723743-144212365 > AG18359.asm.hg38.GYPA_GYPB.bam
samtools index AG18359.asm.hg38.GYPA_GYPB.bam

# AG18359 - permissive mapping
minimap2 -ax asm5 -t 64 --eqx ../references/hg38.noalt.fa /share/dennislab/users/dcsoto/ms_asm/1_asm_shasta/ShastaRun_Jan2023/Assembly.fasta | samtools sort -O BAM - > AG18359.hg38.plotsr.bam
samtools index AG18359.hg38.plotsr.bam
samtools view -Sb AG18359.hg38.plotsr.bam chr4:143723743-144212365 > AG18359.hg38.plotsr.GYPA_GYPB.bam
samtools index AG18359.hg38.plotsr.GYPA_GYPB.bam

# PanTro6 - permissive mapping
minimap2 -ax asm5 -t 64 --eqx ../references/hg38.noalt.fa ../references/panTro6.fa | samtools sort -O BAM - > PanTro6.hg38.plotsr.bam
samtools index PanTro6.hg38.plotsr.bam
samtools view -Sb PanTro6.hg38.plotsr.bam chr4:143723743-144212365 > PanTro6.hg38.plotsr.GYPA_GYPB.bam
samtools index PanTro6.hg38.plotsr.GYPA_GYPB.bam
```

#### 1.1.2 liftOff

Second, we directly lifted over GYP genes to assess their location in each assembly independently of Hg38 structure.

```bash
cd /share/dennislab/users/dcsoto/other/GYPA/annotations
conda activate liftoff

grep "^#\|=GYPA\|=GYPE\|=GYPB" gencode.v43.annotation.gff3 > gencode.v43.annotation.GYP_genes.gff3

# lifting genes from Hg38 to PanTro6
liftoff -copies -sc 0.9 -mm2_options="-t 64 -a --end-bonus 5 --eqx -N 50 -p 0.5 -r 2k -z 5000 -x asm20" -g gencode.v43.annotation.GYP_genes.gff3 -o panTro6.gff ../references/panTro6.fa ../references/hg38.noalt.fa 

liftoff -copies -sc 0.5 -mm2_options="-t 64 -a --end-bonus 5 --eqx -N 50 -p 0.5 -r 2k -z 5000 -x asm20" -g gencode.v43.annotation.GYP_genes.gff3 -o panTro6.2.gff ../references/panTro6.fa ../references/hg38.noalt.fa 

# lifting genes from Hg38 to AG18354
liftoff -copies -sc 0.9 -mm2_options="-t 64 -a --end-bonus 5 --eqx -N 50 -p 0.5 -r 2k -z 5000 -x asm20" -g gencode.v43.annotation.GYP_genes.gff3 -o mPanTro3.hap1.gff ../assemblies/mPanTro3.hap1.20221111.fasta ../references/hg38.noalt.fa
liftoff -copies -sc 0.9 -mm2_options="-t 64 -a --end-bonus 5 --eqx -N 50 -p 0.5 -r 2k -z 5000 -x asm20" -g gencode.v43.annotation.GYP_genes.gff3 -o mPanTro3.hap2.gff ../assemblies/mPanTro3.hap2.20221111.fasta ../references/hg38.noalt.fa

# lifting genes from Hg38 to AG18359
liftoff -copies -sc 0.9 -mm2_options="-t 64 -a --end-bonus 5 --eqx -N 50 -p 0.5 -r 2k -z 5000 -x asm20" -g gencode.v43.annotation.GYP_genes.gff3 -o AG18359.gff /share/dennislab/users/dcsoto/ms_asm/1_asm_shasta/ShastaRun_Jan2023/Assembly.fasta ../references/hg38.noalt.fa

# lifting genes from Hg38 to T2T primates Nov 2023
liftoff -copies -sc 0.9 -mm2_options="-t 64 -a --end-bonus 5 --eqx -N 50 -p 0.5 -r 2k -z 5000 -x asm20" -g gencode.v43.annotation.GYP_genes.gff3 -o mPanTro3.pri.20231122.gff ../assemblies/mPanTro3.pri.cur.20231122.fasta ../references/hg38.noalt.fa
liftoff -copies -sc 0.9 -mm2_options="-t 64 -a --end-bonus 5 --eqx -N 50 -p 0.5 -r 2k -z 5000 -x asm20" -g gencode.v43.annotation.GYP_genes.gff3 -o mPanTro3.alt.20231122.gff ../assemblies/mPanTro3.alt.cur.20231122.fasta ../references/hg38.noalt.fa

liftoff -copies -sc 0.9 -mm2_options="-t 64 -a --end-bonus 5 --eqx -N 50 -p 0.5 -r 2k -z 5000 -x asm20" -g gencode.v43.annotation.GYP_genes.gff3 -o mPanPan1.pri.20231122.gff ../assemblies/mPanPan1.pri.cur.20231122.fasta ../references/hg38.noalt.fa
liftoff -copies -sc 0.9 -mm2_options="-t 64 -a --end-bonus 5 --eqx -N 50 -p 0.5 -r 2k -z 5000 -x asm20" -g gencode.v43.annotation.GYP_genes.gff3 -o mPanPan1.alt.20231122.gff ../assemblies/mPanPan1.alt.cur.20231122.fasta ../references/hg38.noalt.fa
```

### 1.2 Long reads

We directly analyzed long reads to genotype the candidate SNPs.

Obtaining long reads in GYPA/GYPB/GYPE in hg38:
```bash
# AG18359 ONT to Hg38 | Pore chemistry: r941_prom_sup_g5014
cd /share/dennislab/users/dcsoto/other/GYPA/long_reads/AG18359_ont

minimap2 -t 64 -ax map-ont ../references/hg38.noalt.fa /share/dennislab/users/dcsoto/ms_asm/0_ont/AG18359.guppy_5.0.11.fastq.gz | samtools view -F 4 -b - | samtools sort > AG18359.guppy_5.0.11.hg38.srt.bam
samtools index AG18359.guppy_5.0.11.hg38.srt.bam
samtools view -Sb AG18359.guppy_5.0.11.hg38.srt.bam chr4:143723743-144212365 > AG18359.guppy_5.0.11.hg38.srt.GYPA_GYPB.bam
samtools index AG18359.guppy_5.0.11.hg38.srt.GYPA_GYPB.bam

# AG18354 PacBio HiFi to Hg38
cd /share/dennislab/users/dcsoto/other/GYPA/long_reads/AG18354_hifi

wget https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/genomic_data/pacbio_hifi/m54329U_220226_122930.hifi_reads.fastq.gz
wget https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/genomic_data/pacbio_hifi/m54329U_220304_132403.hifi_reads.fastq.gz
wget https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/genomic_data/pacbio_hifi/m64076_210810_005444.hifi_reads.fastq.gz
wget https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/genomic_data/pacbio_hifi/m64076_210813_021703.hifi_reads.fastq.gz
wget https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/genomic_data/pacbio_hifi/m64076_210814_131501.hifi_reads.fastq.gz
wget https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/genomic_data/pacbio_hifi/m64076_210816_001611.hifi_reads.fastq.gz
wget https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/genomic_data/pacbio_hifi/m64076_220303_022914.hifi_reads.fastq.gz

cat m*hifi_reads.fastq.gz > AG18354.hifi_reads.fastq.gz

minimap2 -t 64 -ax map-hifi ../../references/hg38.noalt.fa AG18354.hifi_reads.fastq.gz \
| samtools view -F 4 -b - | samtools sort > AG18354.hifi_reads.hg38.srt.bam
samtools index AG18354.hifi_reads.hg38.srt.bam
samtools view -Sb AG18354.hifi_reads.hg38.srt.bam chr4:143723743-144212365 > AG18354.hifi_reads.hg38.srt.GYPA_GYPB.bam 
samtools index AG18354.hifi_reads.hg38.srt.GYPA_GYPB.bam

# Clint PacBio (Kronenberg et al.)
cd /share/dennislab/users/dcsoto/other/GYPA/long_reads/Clint_ccs

cat download.txt | xargs -n1 -P10 bash -c 'wget $0'

minimap2 -t 64 -ax map-pb ../../references/hg38.noalt.fa Clint_CCS.fastq.gz \
| samtools view -F 4 -b - | samtools sort > Clint_CCS.hg38.srt.bam
#[WARNING] failed to parse the FASTA/FASTQ record next to 'SRR5269278.156764'. Continue anyway.
#[WARNING] failed to parse the FASTA/FASTQ record next to 'SRR5269278.163476'. Continue anyway.
samtools index Clint_CCS.hg38.srt.bam
samtools view -Sb Clint_CCS.hg38.srt.bam chr4:143723743-144212365 > Clint_CCS.hg38.srt.GYPA_GYPB.bam
samtools index Clint_CCS.hg38.srt.GYPA_GYPB.bam

# AG18359 Bionano Genomics
cd /share/dennislab/users/dcsoto/other/GYPA/long_reads/AG18359_bionano

bcftools sort /share/dennislab/users/dcsoto/ms_svs/AG18359/bionano/AG18359.vcf > AG18359_bionano.hg38.vcf
bgzip AG18359_bionano.hg38.vcf; tabix AG18359_bionano.hg38.vcf.gz
bcftools view AG18359_bionano.hg38.vcf.gz -r chr4:143723743-144212365 | bcftools view -f 'PASS,.' | grep -v "IMPRECISE" > AG18359_bionano.hg38.GYPA_GYPB.vcf
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' AG18359_bionano.hg38.GYPA_GYPB.vcf > AG18359_bionano.hg38.GYPA_GYPB.bed
```

We looked for SNPs at candidate SNPs sites as well as paralogous sites in long reads.

**C>A at hg19:chr4:145040845/hg38:chr4:144119692**

We selected a 130 bp region covering GYPA exon 3 for BLAT.
- GYPA: chr4:144119671-144119800 
- GYPE: chr4:143879595-143879721
- GYPB: chr4:144000368-144000494

We aligned DNA from these regions using Mega/MUSCLE to find coordinates of candidate SNPs.
- GYPA hg19:chr4:145040845/hg38:chr4:144119692
- GYPE hg38:chr4:143879616
- GYPB hg38:chr4:144000389

#### 1.2.1 Calling SNVs per contig

We called SNVs and SVs using long reads as contigs in AG18354.

```bash
cd /share/dennislab/users/dcsoto/other/GYPA/long_reads

# selecting reads overlapping candidate allele
bedtools intersect -wa -a AG18354.hifi_reads.hg38.srt.GYPA_GYPB.bam -b candidate_allele.bed \
| samtools view -Sb -F 256 > AG18354.hifi_reads.hg38.srt.GYPA_only.bam

# splitting into individual files
samtools view -H AG18354.hifi_reads.hg38.srt.GYPA_only.bam > AG18354.hifi_reads.hg38.srt.GYPA_only.header.txt
samtools view AG18354.hifi_reads.hg38.srt.GYPA_only.bam > AG18354.hifi_reads.hg38.srt.GYPA_only.sam
split -d -l1 AG18354.hifi_reads.hg38.srt.GYPA_only.sam AG18354.hifi_reads.hg38.srt.GYPA_only_

for num in {00..90} 90{00..59}; do
    cat AG18354.hifi_reads.hg38.srt.GYPA_only.header.txt AG18354.hifi_reads.hg38.srt.GYPA_only_${num} \
    | samtools view -Sb > AG18354.hifi_reads.hg38.srt.GYPA_only_${num}.bam
    samtools index AG18354.hifi_reads.hg38.srt.GYPA_only_${num}.bam
done

# calling SNVs per reads
cd /share/dennislab/users/dcsoto/other/GYPA/long_reads
mkdir -p SNVs_Hg38_AG18354

conda activate assembly
for num in {00..90} 90{00..59}; do
    htsbox pileup -q 0 -evcf ../references/hg38.noalt.fa AG18354_hifi/AG18354.hifi_reads.hg38.srt.GYPA_only_${num}.bam \
    > SNVs_Hg38_AG18354/AG18354.hifi_reads.hg38.srt.GYPA_only.snps_${num}.vcf
    bgzip -f SNVs_Hg38_AG18354/AG18354.hifi_reads.hg38.srt.GYPA_only.snps_${num}.vcf
    tabix -f SNVs_Hg38_AG18354/AG18354.hifi_reads.hg38.srt.GYPA_only.snps_${num}.vcf.gz
done

# merging all SNPs calls into "multi sample" (multi read) VCF files
vcfs=$(ls SNVs_Hg38_AG18354/*.vcf.gz | tr '\n' ' ')
bcftools merge $vcfs -o AG18354.hifi_reads.hg38.srt.GYPA_only.snps.vcf
bgzip AG18354.hifi_reads.hg38.srt.GYPA_only.snps.vcf
tabix AG18354.hifi_reads.hg38.srt.GYPA_only.snps.vcf.gz

# extracting all reads with the candidate alle and  
candidate_reads=$(bcftools view -r chr4:144119692-144119693 AG18354.hifi_reads.hg38.srt.GYPA_only.snps.vcf.gz | bcftools query -f'[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' | grep -w "1/1" | cut -f5 | tr '\n' ' ')

samtools merge -f -o candidate_reads.bam ${candidate_reads}
samtools index candidate_reads.bam
```

```bash
# calling SVs per read
conda activate variants
mkdir -p SVs_Hg38_AG18354

for num in {00..90}; do
    svim-asm haploid --min_sv_size 30 SVs_Hg38_AG18354/${num} AG18354.hifi_reads.hg38.srt.GYPA_only_${num}.bam ../references/hg38.noalt.fa
done

conda activate jasmine
ls SVs_Hg38_AG18354/*/variants.vcf > filelist.txt
jasmine spec_len=10 file_list=filelist.txt out_file=AG18354.hifi_reads.hg38.srt.GYPA_only.svs.vcf
grep -v "SVTYPE=INS" AG18354.hifi_reads.hg38.srt.GYPA_only.svs.vcf > AG18354.hifi_reads.hg38.srt.GYPA_only.svs.dels.vcf

# I manually removed the last variant which is a short-tandem repeat expansion

# To run SVtyper, we neede to add "CIPOS=-100,100;CIEND=-100,100" to each SV entry in vcf file.
awk 'BEGIN{OFS="\t"}{if($0~"^#"){print}else{$8=$8";CIPOS=-100,100;CIEND=-100,100"; print $0}}' AG18354.hifi_reads.hg38.srt.GYPA_only.svs.dels.vcf > AG18354.hifi_reads.hg38.srt.GYPA_only.svs.dels.cipos.vcf
```

> More info about CIPOS in this GitHub issue: https://github.com/hall-lab/svtyper/pull/111

Next, we genotyped variants in `AG18354.hifi_reads.hg38.srt.GYPA_only.svs.dels.cipos.vcf` in Illumina reads (next section).

### 1.3 FastCN

Downloading samples from de Manuel 2016 and Prado-Martinez 2013:
```bash
# Prado-Martinez et al. 2013 samples
cd /share/dennislab/databases/data/great_ape_GP

cat SRR748179_1.fastq.gz SRR748180_1.fastq.gz SRR748181_1.fastq.gz SRR748182_1.fastq.gz SRR748183_1.fastq.gz SRR748184_1.fastq.gz SRR748185_1.fastq.gz SRR748186_1.fastq.gz > Clint_1.fastq.gz
cat SRR748179_2.fastq.gz SRR748180_2.fastq.gz SRR748181_2.fastq.gz SRR748182_2.fastq.gz SRR748183_2.fastq.gz SRR748184_2.fastq.gz SRR748185_2.fastq.gz SRR748186_2.fastq.gz > Clint_2.fastq.gz

cat SRR747951_1.fastq.gz SRR747952_1.fastq.gz SRR747953_1.fastq.gz SRR747954_1.fastq.gz > Bosco_1.fastq.gz
cat SRR747951_2.fastq.gz SRR747952_2.fastq.gz SRR747953_2.fastq.gz SRR747954_2.fastq.gz > Bosco_2.fastq.gz

cat SRR747947_1.fastq.gz SRR747948_1.fastq.gz SRR747949_1.fastq.gz SRR747950_1.fastq.gz > Donald_1.fastq.gz
cat SRR747947_2.fastq.gz SRR747948_2.fastq.gz SRR747949_2.fastq.gz SRR747950_2.fastq.gz > Donald_2.fastq.gz

cat SRR748051_1.fastq.gz SRR748052_1.fastq.gz SRR748053_1.fastq.gz SRR748054_1.fastq.gz SRR748055_1.fastq.gz > Jimmie_1.fastq.gz
cat SRR748051_2.fastq.gz SRR748052_2.fastq.gz SRR748053_2.fastq.gz SRR748054_2.fastq.gz SRR748055_2.fastq.gz > Jimmie_2.fastq.gz

cat SRR748081_1.fastq.gz SRR748082_1.fastq.gz SRR748083_1.fastq.gz SRR748084_1.fastq.gz > Koby_1.fastq.gz
cat SRR748081_2.fastq.gz SRR748082_2.fastq.gz SRR748083_2.fastq.gz SRR748084_2.fastq.gz > Koby_2.fastq.gz

cat SRR748071_1.fastq.gz SRR748072_1.fastq.gz SRR748073_1.fastq.gz SRR748074_1.fastq.gz SRR748075_1.fastq.gz > Clara_1.fastq.gz
cat SRR748071_2.fastq.gz SRR748072_2.fastq.gz SRR748073_2.fastq.gz SRR748074_2.fastq.gz SRR748075_2.fastq.gz > Clara_2.fastq.gz

cat SRR748061_1.fastq.gz SRR748062_1.fastq.gz SRR748063_1.fastq.gz SRR748064_1.fastq.gz SRR748065_1.fastq.gz > Doris_1.fastq.gz
cat SRR748061_2.fastq.gz SRR748062_2.fastq.gz SRR748063_2.fastq.gz SRR748064_2.fastq.gz SRR748065_2.fastq.gz > Doris_2.fastq.gz

cat SRR748066_1.fastq.gz SRR748067_1.fastq.gz SRR748068_1.fastq.gz SRR748069_1.fastq.gz SRR748070_1.fastq.gz > Julie_A959_1.fastq.gz
cat SRR748066_2.fastq.gz SRR748067_2.fastq.gz SRR748068_2.fastq.gz SRR748069_2.fastq.gz SRR748070_2.fastq.gz > Julie_A959_2.fastq.gz

cat SRR748056_1.fastq.gz SRR748057_1.fastq.gz SRR748058_1.fastq.gz SRR748059_1.fastq.gz SRR748060_1.fastq.gz > Vaillant_1.fastq.gz
cat SRR748056_2.fastq.gz SRR748057_2.fastq.gz SRR748058_2.fastq.gz SRR748059_2.fastq.gz SRR748060_2.fastq.gz > Vaillant_2.fastq.gz

cat SRR726275_1.fastq.gz SRR726278_1.fastq.gz SRR726279_1.fastq.gz SRR726280_1.fastq.gz > Andromeda_1.fastq.gz
cat SRR726275_2.fastq.gz SRR726278_2.fastq.gz SRR726279_2.fastq.gz SRR726280_2.fastq.gz > Andromeda_2.fastq.gz

cat SRR726352_1.fastq.gz SRR726353_1.fastq.gz SRR726354_1.fastq.gz SRR726355_1.fastq.gz SRR726356_1.fastq.gz SRR726357_1.fastq.gz SRR726358_1.fastq.gz SRR726359_1.fastq.gz > Bwambale_1.fastq.gz
cat SRR726352_2.fastq.gz SRR726353_2.fastq.gz SRR726354_2.fastq.gz SRR726355_2.fastq.gz SRR726356_2.fastq.gz SRR726357_2.fastq.gz SRR726358_2.fastq.gz SRR726359_2.fastq.gz > Bwambale_2.fastq.gz

cat SRR747955_1.fastq.gz SRR747956_1.fastq.gz SRR747957_1.fastq.gz SRR747958_1.fastq.gz > Harriet_1.fastq.gz
cat SRR747955_2.fastq.gz SRR747956_2.fastq.gz SRR747957_2.fastq.gz SRR747958_2.fastq.gz > Harriet_2.fastq.gz

cat SRR726360_1.fastq.gz SRR726381_1.fastq.gz SRR726395_1.fastq.gz SRR726402_1.fastq.gz SRR726403_1.fastq.gz SRR726404_1.fastq.gz SRR726405_1.fastq.gz SRR726406_1.fastq.gz SRR726407_1.fastq.gz > Kidongo_1.fastq.gz
cat SRR726360_2.fastq.gz SRR726381_2.fastq.gz SRR726395_2.fastq.gz SRR726402_2.fastq.gz SRR726403_2.fastq.gz SRR726404_2.fastq.gz SRR726405_2.fastq.gz SRR726406_2.fastq.gz SRR726407_2.fastq.gz > Kidongo_2.fastq.gz

cat SRR726408_1.fastq.gz SRR726409_1.fastq.gz SRR726410_1.fastq.gz SRR726411_1.fastq.gz SRR726412_1.fastq.gz SRR726413_1.fastq.gz SRR726415_1.fastq.gz SRR726416_1.fastq.gz > Nakuu_1.fastq.gz
cat SRR726408_2.fastq.gz SRR726409_2.fastq.gz SRR726410_2.fastq.gz SRR726411_2.fastq.gz SRR726412_2.fastq.gz SRR726413_2.fastq.gz SRR726415_2.fastq.gz SRR726416_2.fastq.gz > Nakuu_2.fastq.gz

cat SRR726233_1.fastq.gz SRR726241_1.fastq.gz SRR726242_1.fastq.gz SRR726243_1.fastq.gz > Vincent_1.fastq.gz
cat SRR726233_2.fastq.gz SRR726241_2.fastq.gz SRR726242_2.fastq.gz SRR726243_2.fastq.gz > Vincent_2.fastq.gz

cat SRR748121_1.fastq.gz SRR748122_1.fastq.gz SRR748123_1.fastq.gz SRR748124_1.fastq.gz > Akwaya-Jean_1.fastq.gz
cat SRR748121_2.fastq.gz SRR748122_2.fastq.gz SRR748123_2.fastq.gz SRR748124_2.fastq.gz > Akwaya-Jean_2.fastq.gz

cat SRR748125_1.fastq.gz SRR748126_1.fastq.gz SRR748127_1.fastq.gz SRR748128_1.fastq.gz > Banyo_1.fastq.gz
cat SRR748125_2.fastq.gz SRR748126_2.fastq.gz SRR748127_2.fastq.gz SRR748128_2.fastq.gz > Banyo_2.fastq.gz

cat SRR748129_1.fastq.gz SRR748130_1.fastq.gz SRR748131_1.fastq.gz SRR748132_1.fastq.gz > Basho_1.fastq.gz
cat SRR748129_2.fastq.gz SRR748130_2.fastq.gz SRR748131_2.fastq.gz SRR748132_2.fastq.gz > Basho_2.fastq.gz

cat SRR748133_1.fastq.gz SRR748134_1.fastq.gz SRR748135_1.fastq.gz SRR748136_1.fastq.gz SRR748137_1.fastq.gz > Damian_1.fastq.gz
cat SRR748133_2.fastq.gz SRR748134_2.fastq.gz SRR748135_2.fastq.gz SRR748136_2.fastq.gz SRR748137_2.fastq.gz > Damian_2.fastq.gz

cat SRR748138_1.fastq.gz SRR748139_1.fastq.gz SRR748140_1.fastq.gz SRR748141_1.fastq.gz SRR748142_1.fastq.gz > Julie_LWC21_1.fastq.gz
cat SRR748138_2.fastq.gz SRR748139_2.fastq.gz SRR748140_2.fastq.gz SRR748141_2.fastq.gz SRR748142_2.fastq.gz > Julie_LWC21_2.fastq.gz

cat SRR748143_1.fastq.gz SRR748144_1.fastq.gz SRR748145_1.fastq.gz SRR748146_1.fastq.gz > Kopongo_1.fastq.gz
cat SRR748143_2.fastq.gz SRR748144_2.fastq.gz SRR748145_2.fastq.gz SRR748146_2.fastq.gz > Kopongo_2.fastq.gz

cat SRR748147_1.fastq.gz SRR748148_1.fastq.gz SRR748149_1.fastq.gz SRR748150_1.fastq.gz SRR748151_1.fastq.gz > Koto_1.fastq.gz
cat SRR748147_2.fastq.gz SRR748148_2.fastq.gz SRR748149_2.fastq.gz SRR748150_2.fastq.gz SRR748151_2.fastq.gz > Koto_2.fastq.gz

cat SRR748152_1.fastq.gz SRR748153_1.fastq.gz SRR748154_1.fastq.gz SRR748155_1.fastq.gz > Paquita_1.fastq.gz
cat SRR748152_2.fastq.gz SRR748153_2.fastq.gz SRR748154_2.fastq.gz SRR748155_2.fastq.gz > Paquita_2.fastq.gz

cat SRR748156_1.fastq.gz SRR748157_1.fastq.gz SRR748158_1.fastq.gz SRR748159_1.fastq.gz SRR748160_1.fastq.gz > Taweh_1.fastq.gz
cat SRR748156_2.fastq.gz SRR748157_2.fastq.gz SRR748158_2.fastq.gz SRR748159_2.fastq.gz SRR748160_2.fastq.gz > Taweh_2.fastq.gz

cat SRR748161_1.fastq.gz SRR748162_1.fastq.gz SRR748163_1.fastq.gz SRR748164_1.fastq.gz > Tobi_1.fastq.gz
cat SRR748161_2.fastq.gz SRR748162_2.fastq.gz SRR748163_2.fastq.gz SRR748164_2.fastq.gz > Tobi_2.fastq.gz

# Downloading additional samples from ENA
cd /share/dennislab/users/dcsoto/other/GYPA/download
/share/dennislab/users/dcsoto/Miniconda3/envs/smk/bin/snakemake --snakefile download.smk --config input=samples.tsv2 -p -j 5

# Downloading additional samples from Harrison/Bioshare and GDrive
cd /share/dennislab/users/dcsoto/other/GYPA/crams

conda activate assembly # newest version of samtools required
/share/dennislab/users/dcsoto/Miniconda3/envs/smk/bin/snakemake --snakefile snake-cram2fastq.smk --config input=input.tsv -p -j 64 -n
```

All samples from de Manuel 2016 and Prado-Martinez 2013 were added to  `/share/dennislab/databases/data/great_ape_GP/` and linked files to `fastqs` folder.
```bash
cd /share/dennislab/users/dcsoto/other/GYPA/fastqs

# AG18359
ln -s /share/dennislab-backedup/illumina/wgs/chimp/AG18359/hwftp.novogene.com/C202SC19030370.DHG.PE150.20190412.P202SC19030370-01-01/data_release/C202SC19030370/raw_data/AG18359/AG18359_USD16090428L_HKGF3DSXX_L4_1.fq.gz AG18359_1.fastq.gz
ln -s /share/dennislab-backedup/illumina/wgs/chimp/AG18359/hwftp.novogene.com/C202SC19030370.DHG.PE150.20190412.P202SC19030370-01-01/data_release/C202SC19030370/raw_data/AG18359/AG18359_USD16090428L_HKGF3DSXX_L4_2.fq.gz AG18359_2.fastq.gz

# S003641
ln -s /share/dennislab-backedup/illumina/wgs/chimp/S003641/S003641_S11_L002_R1_001.fastq.gz S003641_1.fastq.gz
ln -s /share/dennislab-backedup/illumina/wgs/chimp/S003641/S003641_S11_L002_R2_001.fastq.gz S003641_2.fastq.gz

# de Manuel 2016 and Prado-Martinez 2013 samples
ln -s /share/dennislab/databases/data/great_ape_GP/*gz .
```

We genotyped the copy number of GYPA full lenght transcript and the copy number of the 16 kbp deletion and 3 kbp deletion regions. Considering that QuickMer2 is noisy in this region (and PanTro6 is not a good reference), we used gene-family copy numbers (FastCN/WSSD).

We re-run FastCN in samples from the Great Ape Project:
```bash
cd /share/dennislab/users/dcsoto/other/GYPA/fastCN/

conda activate snakecn
export PATH="/share/dennislab/programs/fastCN:/share/dennislab/programs/mrsfast/:$PATH"

/share/dennislab/users/dcsoto/Miniconda3/envs/smk/bin/snakemake --snakefile snake-fastcn.smk --config reference_path=GRCh38_BSM_WMDUST windows_size=1kb input=input.tsv data=../fastqs -p -j 32
```

We converted deletions into a bed file to genotype them and added additional regions of interest (GRCh38):
```bash
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\n' ../long_reads/AG18354.hifi_reads.hg38.srt.GYPA_only.svs.dels.cipos.vcf > deletions.bed
cat deletions.bed <(echo -e "chr4\t144115772\t144140695\tGYPA_dup") <(echo -e "chr4\t144109298\t144115775\tGYPA_unique") > regions.bed
```

> GYPA stringent unique region: hg38: chr4:144101660-144111660

We genotyped CN of these deletions in these individuals:
```bash
conda activate python3
python3 genotype_cn_parallel.py --path results/windows --genes regions.bed --output regions.cn.tsv -t 64
python3 genotype_cn_parallel.py --path results/windows --genes segdups.bed --output segdups.cn.tsv -t 64
```

BigBed files were stored in BioShare: `https://bioshare.bioinformatics.ucdavis.edu/bioshare/view/Dennis_UCSC_Browser_Files/dcsoto/GYPA/FastCN/`

Files were uplaoded to the UCSC Genome Browser:
```bash
for file in *CN.bb; do name=$(basename -s .depth.1kb.bed.CN.bb $file); echo "track type=bigBed name=\"$name\" description=\"FastCN CN ${name}\" itemRgb=On bigDataUrl=\"https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/cpqqdfge5lfvovq/dcsoto/GYPA/FastCN/${name}.depth.1kb.bed.CN.bb\""; done
```

## 2. PanTro6 as reference

### 2.1 Long reads

```bash
cd /share/dennislab/users/dcsoto/other/GYPA/GYPA_long

# AG18359 ONT to PanTro6 | Pore chemistry: r941_prom_sup_g5014
minimap2 -t 64 -ax map-ont ../references/panTro6.fa /share/dennislab/users/dcsoto/ms_asm/0_ont/AG18359.guppy_5.0.11.fastq.gz | samtools view -F 4 -b - | samtools sort > AG18359.guppy_5.0.11.panTro6.srt.bam
samtools index AG18359.guppy_5.0.11.panTro6.srt.bam
samtools view -Sb AG18359.guppy_5.0.11.panTro6.srt.bam chr4:141934474-142246523 > AG18359.guppy_5.0.11.panTro6.srt.GYPA_GYPB.bam
samtools index AG18359.guppy_5.0.11.panTro6.srt.GYPA_GYPB.bam

# AG18354 PacBio HiFi to PanTro6
minimap2 -t 64 -ax map-ont ../references/panTro6.fa AG18354.hifi_reads.fastq.gz | samtools view -F 4 -b - | samtools sort > AG18354.hifi_reads.panTro6.srt.bam
samtools index AG18354.hifi_reads.panTro6.srt.bam
samtools view -Sb AG18354.hifi_reads.panTro6.srt.bam chr4:141934474-142246523 > AG18354.hifi_reads.panTro6.srt.GYPA_GYPB.bam
samtools index AG18354.hifi_reads.panTro6.srt.GYPA_GYPB.bam
```

### 2.2 QuickMer2

First, we obtained SegDups in PanTro6:
```bash
module load samtools/1.9
samtools faidx /share/dennislab/databases/assemblies/panTro6/panTro6.fa

cd /share/dennislab/programs/sedef
./sedef.sh -o /share/dennislab/databases/assemblies/panTro6/sedef -j5 /share/dennislab/databases/assemblies/panTro6/panTro6.fa
```

SDs available at `/share/dennislab/databases/assemblies/panTro6/sedef/final.bed`

We converted SDs from BEDPE to BED:
```bash
cd /share/dennislab/databases/assemblies/panTro6/sedef/
#python3 /share/dennislab/programs/sedef/sedef_to_bed_mod.py final.bed panTro6.SDs.bed panTro6.SDs.lowid.bed
python3 sedef_to_bed_mod.py final.bed panTro6.SDs.bed panTro6.SDs.lowid.bed # we removed first line in panTro6.SDs.bed
```

> Additionally, we run BISER as follow: (but we didn't use this output)
> ```bash
> conda create -n biser
> conda activate biser
> conda install python=3.7
> pip install biser
> conda install -c bioconda samtools
> cd /share/dennislab/databases/assemblies/panTro6/biser
> conda activate biser 
> biser --output out -t 10 --temp temp /share/dennislab/databases/assemblies/panTro6/panTro6.fa 
> ```

Then, we prepared control regions (regions known to be diploid copy-number fixed, i.e., removing sex chromosomes, mithocondrial DNA, and segdups):
```bash
cd /share/dennislab/databases/data/quickmer-chimp

bedtools complement -i /share/dennislab/databases/assemblies/panTro6/sedef/panTro6.SDs.bed -g /share/dennislab/databases/assemblies/panTro6/panTro6.fa.fai | grep -v "chrM\|chrX\|chrY" > control/panTro6.control.bed
```

Finally, we prepared samples:
```bash
cd /share/dennislab/databases/data/quickmer-chimp/fastqs

cat /share/dennislab-backedup/illumina/wgs/chimp/AG18359/hwftp.novogene.com/C202SC19030370.DHG.PE150.20190412.P202SC19030370-01-01/data_release/C202SC19030370/raw_data/AG18359/AG18359_USD16090428L_HKGF3DSXX_L4_1.fq.gz /share/dennislab-backedup/illumina/wgs/chimp/AG18359/hwftp.novogene.com/C202SC19030370.DHG.PE150.20190412.P202SC19030370-01-01/data_release/C202SC19030370/raw_data/AG18359/AG18359_USD16090428L_HKGF3DSXX_L4_2.fq.gz > AG18359.fastq.gz

cat /share/dennislab-backedup/illumina/wgs/chimp/S003641/S003641_S11_L002_R1_001.fastq.gz /share/dennislab-backedup/illumina/wgs/chimp/S003641/S003641_S11_L002_R2_001.fastq.gz > S003641.fastq.gz
```

Running QuicKmer2 in chimpanzees:
```bash
cd /share/dennislab/databases/data/quickmer-chimp

conda activate variants # python3 with pandas and argparse, and samtools + bedToBigBed
export PATH=/share/dennislab/databases/data/quickmer-t2t/QuicK-mer2/:$PATH

/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake  \
--snakefile quickmer2.smk \
--config reference=reference/panTro6.fa control=control/panTro6.control.bed \
input=fastqs outdir=results_w500 distance=1 window=500 \
-p -j 20
```
