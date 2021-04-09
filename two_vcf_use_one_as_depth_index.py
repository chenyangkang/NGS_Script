#!/usr/bin/python
import os
import pandas as pd
import numpy as np
import sys
import re
import shutil

def open_vcf1(input):
    f=open(input,"r")
    global vcf1
    vcf1=f.read()
    vcf1=vcf1.strip().split("\n")
    f.close()

def open_vcf2(input):
    f=open(input,"r")
    global vcf2
    vcf2=f.read()
    vcf2=vcf2.strip().split("\n")
    f.close()

def get_qualified_SNPs(vcf1_input,vcf2_input,depth,output):
    open_vcf1(vcf1_input)
    open_vcf2(vcf2_input)
    output=open("%s"%(output),"w")
    dict={}
    line_count=0
    for line in vcf1:
        if line.startswith("#CHROM"):
            sample_num=len(line.split("\t"))-9
        else:
            if not line.startswith("#"):
                line1=line.split("\t")
                DP_pos=line1[8].split(":").index("DP")
                flag=1
                for sample in range(1,sample_num+1):
#                    print("before  %s"%(line1[8+sample].split(":")[DP_pos]))
                    if int(line1[8+sample].split(":")[DP_pos])>= int(depth):
#                        print("after     %s"%(line1[8+sample].split(":")[DP_pos]))
                        continue
                    else:
                        flag=flag-2
                if flag>=1:
                    line_count+=1
                    try: 
                        dict[line1[0]].append(line1[1])
                    except:
                        dict[line1[0]]=[]
                        dict[line1[0]].append(line1[1])

    for line in vcf2:
        if line.startswith("#"):
            output.write("%s\n"%(line))
        else:
            line_parse=line.split("\t")
            if line_parse[1] in dict[line_parse[0]]:
                output.write("%s\n"%(line))
    
    output.close()

if __name__ == "__main__":
    if len(sys.argv) == 5:
        get_qualified_SNPs(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    else:
        print("Usage: get_qualified_SNPs.py vcf1_input(as_the_larger_one_with_depth_info) vcf2_input depth output")
        sys.exit(0)
        
        
