#!/usr/bin/env python3

# Import libraries
import pandas as pd
import sys
import argparse



def join_fasta(args):
    file_in_prrt = open(args.prrt)
    file_in_int = open(args.int)
    file_out = open(args.output, "w")
    scount = None
    before_prrt_region = ""
    prrt_region = ""
    between_prrt_int_region = ""
    int_region = ""
    after_int_region = ""
    ref = ""
    for position in range(1, 9719):
        if position < 2277:
             before_prrt_region =  before_prrt_region + "N"
        elif position >= 2277 and position < 3507:
             for line in file_in_prrt:
                if ">" in line:
                    scount = line.split("_")[0]
                else:
                    for nucleotide in line:
                        prrt_region = prrt_region + nucleotide
                        prrt_region = prrt_region.strip("\n")
        elif position > 3506 and position < 4230:
            between_prrt_int_region =  between_prrt_int_region + "N"
        elif position > 4229 and position < 5064:
            for line in file_in_int:
                if line[0] != ">":
                    for nucleotide in line:
                        int_region = int_region + nucleotide
                        int_region = int_region.strip("\n")
        else:
            after_int_region =  after_int_region + "N"
            


    
                 
    file_in_prrt.close()
    file_in_int.close()
    ref = before_prrt_region + prrt_region + between_prrt_int_region + int_region + after_int_region
    print(len(ref)) 
    print(len(prrt_region))
    print(len(int_region))
    file_out.write(scount + "_prrt" + "_int")
    file_out.write("\n")
    file_out.write(ref)
    file_out.write("\n")
    file_out.close()
 


def main():
    parser=argparse.ArgumentParser(description = "Merge PRRT and INT relative to position in HXB2")
    parser.add_argument("-p", help="PRRT fasta input file", dest = "prrt", type = str, required=True)
    parser.add_argument("-i", help="INT fasta input file", dest="int", type = str, required=True)
    parser.add_argument("-o", help ="fastq output filename", dest = "output", type = str, required=True)
    parser.set_defaults(func = join_fasta)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
	main()