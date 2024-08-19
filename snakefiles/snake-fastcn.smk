"""
Snakefile FastCN-lite
"""

import os
import re
import pandas as pd


# CONFIG FILE #


REFERENCE_PATH = config['reference_path']
WINDOWSIZE = config['windows_size']
INPUT = config['input']
DATA = config['data']

input_file = pd.read_csv(INPUT, delimiter = '\t', names = ['sample','forward','reverse'])
samples = list(set(input_file['sample']))
samples_dict = input_file.set_index('sample').to_dict('index')


#### WORKFLOW TARGET ####


rule all:
  input:
    expand("results/bigBed/{smp}.depth."+WINDOWSIZE+".bed.CN.bb", smp = samples)


#### REFERENCE PROCESSING ####


rule gc_correction_ref:
  input:
    fasta = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST_unmasked.fa",
    exclude = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST.exclude.bed.sort2",
    gaps = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST.gaps.bed.slop36.sorted.merged.sort2"
  output:
    bin = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST.GC_control.bin"
  shell:
    '''
    GC_control_gen {input.fasta} {input.exclude} {input.gaps} 400 {output.bin}
    '''


rule index_masked:
  input:
    fasta = REFERENCE_PATH+"/ref-WMDUST-masked/GRCh38_BSM_WMDUST_masked.fa"
  output:
    index = REFERENCE_PATH+"/ref-WMDUST-masked/GRCh38_BSM_WMDUST_masked.fa.index"
  shell:
    '''
    mrsfast --index {input.fasta}
    '''
    

#### MAPPING TO REFERENCE ####


def get_forward(wildcards):
  return DATA+"/"+samples_dict.get(wildcards.smp)['forward']
  

def get_reverse(wildcards):
  return DATA+"/"+samples_dict.get(wildcards.smp)['reverse']


rule mapping:
  input:
    r1 = get_forward,
    r2 = get_reverse,
    ref = REFERENCE_PATH+"/ref-WMDUST-masked/GRCh38_BSM_WMDUST_masked.fa",
    index = REFERENCE_PATH+"/ref-WMDUST-masked/GRCh38_BSM_WMDUST_masked.fa.index"
  output:
    temp("results/mapping/{smp}.sam.gz")
  params:
    "results/mapping/{smp}"
  threads: 16
  shell:
    '''
    extract-from-fastq36-pair.py --in1 {input.r1} --in2 {input.r2} | mrsfast --search {input.ref} --seq /dev/fd/0 --disable-nohits --mem 16 --threads {threads} -e 2 --outcomp -o {params}
    '''


#### GC-CORRECTED READ-DEPTH ####


rule gc_correction_sam:
  input:
    alignment = "results/mapping/{smp}.sam.gz",
    reference = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST_unmasked.fa.fai",
    gccontrol = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST.GC_control.bin"
  output:
    bpdepth = temp("results/binary/{smp}.bin.gz")
  params:
    "results/binary/{smp}"
  shell:
    '''
    zcat {input.alignment} | SAM_GC_correction {input.reference} {input.gccontrol} /dev/fd/0 {params}
    gzip {params}.bin
    '''


#### CONVERT READ-DEPTH FROM BP TO WINDOW ####


rule bp_to_windows:
  input:
    bpdepth = "results/binary/{smp}.bin.gz",
    chromlen = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST_unmasked.fa.fai", 
    windows = REFERENCE_PATH+"/windows-WMDUST/GRCh38_BSM_WMDUST."+WINDOWSIZE+".bed"
  output:
    windepth = "results/windows/{smp}.depth."+WINDOWSIZE+".bed"
  shell:
    '''
    perbp-to-windows.py --depth {input.bpdepth} --out {output.windepth} --chromlen {input.chromlen} --windows {input.windows}
    '''


#### READ-DEPTH TO COPY NUMBER ####


rule depth_to_cn:
  input:
     depth = "results/windows/{smp}.depth."+WINDOWSIZE+".bed",
     auto = REFERENCE_PATH+"/windows-WMDUST/GRCh38_BSM_WMDUST."+WINDOWSIZE+".autoControl.bed",
     chrx = REFERENCE_PATH+"/windows-WMDUST/GRCh38_BSM_WMDUST."+WINDOWSIZE+".chrXnonParControl.bed"
  output:
     "results/windows/{smp}.depth."+WINDOWSIZE+".bed.CN.bed"
  shell:
    '''
    depth-to-cn.py --in {input.depth} --autocontrol {input.auto} --chrX {input.chrx}
    '''


#### BIGBED CONVERSION ####


rule bed2bigBed:
  input:
    bedGraph = "results/windows/{smp}.depth."+WINDOWSIZE+".bed.CN.bed",
    chromsizes = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM.chromsizes" 
  output:
    bed9 = temp("results/bigBed/{smp}.depth."+WINDOWSIZE+".bed.CN.bed9"),
    sorted = temp("results/bigBed/{smp}.depth."+WINDOWSIZE+".bed.CN.srt.bed9"),
    bigbed = "results/bigBed/{smp}.depth."+WINDOWSIZE+".bed.CN.bb"
  priority: 100
  shell:
    '''
    python scripts/bedToBed9.py {input.bedGraph} {output.bed9} 
    sort -k1,1 -k2,2n {output.bed9} > {output.sorted}
    bedToBigBed -type=bed9 {output.sorted} {input.chromsizes} {output.bigbed}  
    '''


