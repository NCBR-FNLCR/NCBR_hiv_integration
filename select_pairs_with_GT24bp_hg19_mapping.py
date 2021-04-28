#!/usr/env/bin python
import sys
import re
import os

#start with teh filtered R2 sam generated with the pairwise2_v0.6_KiteLenti_specific.py script
#make sam out of unfiltered R1 sam to use as input to this script

def get_R2_read_names(R2_sam):
    R2_read_names=[]
    with open(R2_sam) as R2_sam:
        for line in R2_sam.readlines()[29:]:
            read_name=line.split('\t')[0]
            R2_read_names.append(read_name)
    return set(R2_read_names)
    
def are_GT20bps_matching_in_R1(R2_read_names,R1_sam):
    filtered_R1_read_names=[]
    with open(R1_sam) as R1_sam:
        for line in R1_sam.readlines()[29:]:
            read_name=line.split('\t')[0]
            if not read_name in R2_read_names:
                continue
            CIGAR=re.split('\t| ',line)[5]
            if len(CIGAR.split('M'))>2:
                # print "This cigar string has >1 matching sections:"
                # print read_name
                for section in CIGAR.split('M'):
                    # print section
                    try:
                        length_matching=int(section)
                        if length_matching>24:
                            filtered_R1_read_names.append(read_name)
                            break
                    except:
                        pass
                        # print "this section didn't have >24 bases matching:",section
            else:
                matching=CIGAR.split('M')[0].split('S')[-1].split('H')[-1]
            try:
                if int(matching)>24:
                    # print "this section of an R1 matches >24 bases"
                    # print read_name
                    # print matching
                    filtered_R1_read_names.append(read_name)
            except:
                pass
                # print "this cigar string couldn't be parsed"
                # print read_name
                # print CIGAR
                
    return set(filtered_R1_read_names)
    
def create_further_filtered_R2_sam(filtered_R1_read_names,R2_sam):
    sample_name=os.path.basename(R2_sam).split('_')[0]
    outfile=open(sample_name+'_further_filtered_R2.sam','w')
    with open(R2_sam) as R2_sam:
        for line in R2_sam.readlines():
            # print line
            if line.startswith('@'):
                outfile.write(line)
                continue
            read_name=line.split('\t')[0]
            if read_name in filtered_R1_read_names:
                outfile.write(line)
    outfile.close()

def main(args):
    
    print """Enter the following inputs:
    1) unfiltered R1 sam
    2) filtered R2 sam"""

    R2_read_names=get_R2_read_names(sys.argv[2])
    filtered_R1_read_names=are_GT20bps_matching_in_R1(R2_read_names,sys.argv[1])
    create_further_filtered_R2_sam(filtered_R1_read_names,sys.argv[2])
    

if __name__ == "__main__":
    main(sys.argv)