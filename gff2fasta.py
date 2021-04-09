#!/usr/bin/python
import os
import pandas as pd
import numpy as np
import sys
import re
import shutil

def open_gff(input):
    f=open(input,"r")
    global gff_all
    gff_all=f.read()
    f.close()

def open_fna(input2):
    f=open(input2,"r")
    global fna_all
    fna_all=f.read()
    f.close()

def get_regions(target_region,gff,fna):
    cwd=os.getcwd()
    if os.path.exists(cwd+'/gene'):
        shutil.rmtree(cwd+'/gene')
    os.mkdir(cwd+'/gene')
    open_gff(gff)
    open_fna(fna)
    dict={}
    line_count=0
    gff_parse=gff_all.strip().split("\n")
    for line in gff_parse:
        if not line.startswith("#"):
            each_line=line.strip().split("\t")
            if len(each_line)!=0:
                if each_line[2]==str(target_region):
                    gene_anno_line=each_line[8].split(";")
                    if "ID=gene" in gene_anno_line[0]:
                        line_count+=1
                        dict[line_count]=[]
                        dict[line_count].append(each_line[0])
                        dict[line_count].append(each_line[3])
                        dict[line_count].append(each_line[4])
                        gene_name=gene_anno_line[1].split("=")[1]
                        if "/" in gene_name:
                            gene_name.replace("/","&")
                        dict[line_count].append(gene_name)
        else:
            continue

    fna_parse=fna_all.strip().split(">")[1:]
    for chr in fna_parse:
        chr_line=chr.split("\n")
        first_line=chr_line[0].split(" ")
        chr_name=first_line[0]
#        print(first_line)
        species_name="%s_%s"%(first_line[1],first_line[2])
        for gene in dict.keys():
            if str(dict[gene][0])==str(chr_name):
                chr_line=chr_line[1:]
                chr_one_line=chr.replace("\n","")
                gene_sequence="\n".join(re.findall(r'.{1,80}',chr_one_line[int(dict[gene][1])-1:int(dict[gene][2])]))
                output=open("./gene/%s.gene"%(dict[gene][3]),"w")
                output.write(">%s\n"%(str(dict[gene][3])))
                output.write("%s\n"%(gene_sequence))
                output.close()

if __name__ == "__main__":
    if len(sys.argv) == 4:
        get_regions(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print("Usage: get_region.py target_regions gff_file genome_file(fna) ('CDS' 'gene' 'mrna' 'exon')")
        sys.exit(0)


### Copy right: This script is written by Chen Yangkang, master student from Institute of Zoology, Chinese Academy of Science.
#  Please donnot remove this block.
# 
# ###
