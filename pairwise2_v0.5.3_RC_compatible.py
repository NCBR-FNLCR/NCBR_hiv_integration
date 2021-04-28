#!//anaconda/bin/python

#Want to determine whether a sequencing read is a mispriming event vs. an actual HIV integration event in the human genome
#After obtaining R2 reads that contain one of 16 possible HIV GSP2 seqs (ie filter out control reads), align said reads to HIV consensus seqs
#Use an alignment scores between HIV consensus seqs and R2s to filter out possible mispriming events as opposed to actual integration events

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
from Bio.Seq import Seq
import subprocess

CONSENSUS_A='TGGATGGGTTAATTtACTCcaaGAaAAGACAAGAAATCCTTGATCTGTGGGTCTAtaACACACAAGGaTtCTTCCCtGATTGGCAGAAtTACACACCAGGGCCAGGgAccAGATtCCC'
CONSENSUS_B='TGGAaGGGcTAaTttacTCcCaaaaaAGaCAaGAtATcCTTGATcTgTGGatctAccACACACAAGGCTacTTCCCtGAtTgGCAGAAcTACACAcCAGGgccagGgatcAGaTatCC'
CONSENSUS_D='TGGAAGGGCTAATTTGGTC?AAaAa?AGACAAGAgATCCTTGAtCTTTGGGTCTACAACACACAAGGCATCTTCCCtGATTGGCAgAACTACACACCAGGGCCAGGGATCAGATAtCC'
CONSENSUS_O='TGGA?GGGTTAATTTACTCCCATAA?AGAGCAGAAATCCTGGATCT?TGG?T?TAT?ACACTCAGGGATTCTTCCCTGATTGGCAG???TACACACC?GGACCAGGA?C?AG?TTCCC'
CONSENSUS_U='TGGATGGG?TA?TTT??TCC???AAAAGACAAGAAATCCTTGA?CTGTGGGTCTA?CA?ACACAAGGCT?CTTCCCTGATTGGCA?AA?TACACACCAGGGCCAGGGA??AGAT?CCC'
CONSENSUS_CPZ='TGGAAGGGTTAGTTTACTCCAGGAGAAGACAAGAGATCCTTGACCTCTGGGTCTATCACACACAAGGCTTCTTCCCTGACTGGCAGAACTACACAACAGGACCAGGAACAAGATTCCC'

Cons=[CONSENSUS_A,CONSENSUS_B,CONSENSUS_D,CONSENSUS_O,CONSENSUS_U,CONSENSUS_CPZ]

Con_names=['ConA','ConB','ConD','ConO','ConO','CCPZ']

GSP2='tcwgggaagtakccttgwgtrt'.upper()
GSP2_seq=Seq(GSP2)
GSP2_revcomp=str(GSP2_seq.reverse_complement())
all_GSP2s=['TCAGGGAAGTAGCCTTGAGTAT','TCAGGGAAGTAGCCTTGAGTGT','TCAGGGAAGTAGCCTTGTGTAT','TCAGGGAAGTAGCCTTGTGTGT','TCAGGGAAGTATCCTTGAGTAT','TCAGGGAAGTATCCTTGAGTGT','TCAGGGAAGTATCCTTGTGTAT','TCAGGGAAGTATCCTTGTGTGT','TCTGGGAAGTAGCCTTGAGTAT','TCTGGGAAGTAGCCTTGAGTGT','TCTGGGAAGTAGCCTTGTGTAT','TCTGGGAAGTAGCCTTGTGTGT','TCTGGGAAGTATCCTTGAGTAT','TCTGGGAAGTATCCTTGAGTGT','TCTGGGAAGTATCCTTGTGTAT','TCTGGGAAGTATCCTTGTGTGT']

Rev_Comps_GSP2s=[]

for gsp2 in all_GSP2s:
    gsp2=gsp2.upper()
    seq=Seq(gsp2)
    revcomp=str(seq.reverse_complement())
    Rev_Comps_GSP2s.append(revcomp)

Rev_Comps_Consensi=[]

for Con in Cons:
    Con=Con.upper()
    seq=Seq(Con)
    Rev_Comps_Consensi.append(seq.reverse_complement())

#samtools view -h -o out.sam in.bam    
cmd1=["samtools", "view", "-h", "-o",sys.argv[1]+".sam",sys.argv[1]]
c1=subprocess.Popen(cmd1)
c1.communicate()


with open(sys.argv[1]+".sam",'r') as sam:
    outfile_name=sys.argv[1].split('.')[0]+'HIV_aligned.sam'
    outfile=open(outfile_name,'w')
    for line in sam:
        if not line.startswith('@'):
            Read2=line.split()[9]
            if any(substring in Read2 for substring in Rev_Comps_GSP2s):
                strand='-'
                Consensus_Seqs_correct_orientation=Cons
                GSP2=GSP2_revcomp
            elif any(substring in Read2 for substring in all_GSP2s):
                strand='+'
                Consensus_Seqs_correct_orientation=Rev_Comps_Consensi
                GSP2='tcwgggaagtakccttgwgtrt'.upper()
            else:
                continue
        else:
            outfile.write(line)
            continue


        best_score=0
        all_alignments=[]
        all_matches=[]
        

        for consensus in Consensus_Seqs_correct_orientation:
            alignments = pairwise2.align.localxs(consensus,Read2,-1,-1)
            
            # print alignments
            gsp2_alignments=pairwise2.align.localxs(GSP2,Read2,-1,-1)

            #now get the best alignments only
#            best_score=0
            best_gsp2_score=0

            for a in alignments:
                al1,al2,score,start,end=a
                if score>best_score:
                    best_score=score
                    best_alignments=a
                    best_start=start-1
        
            for g in gsp2_alignments:
                al1,al2,gsp2_score,start,end=g
                if gsp2_score>best_gsp2_score:
                    best_gsp2_score=gsp2_score
                    best_gsp2_alignments=g
                    best_gsp2_start=start-1

            match = []
            match_gsp2=[]

            for a, b in zip(best_alignments[0],best_alignments[1]):
                if a == b:
                    match.append('|')
                else:
                    match.append(' ')
            
            for c, d in zip(gsp2_alignments[0][0],gsp2_alignments[0][1]):
                if c == d:
                    match_gsp2.append('|')
                elif c=='W' and d =='A':
                    match_gsp2.append('|')
                elif c=='W' and d =='T':
                    match_gsp2.append('|')
                elif c=='K' and d =='G':
                    match_gsp2.append('|')
                elif c=='K' and d =='T':
                    match_gsp2.append('|')
                elif c=='R' and d =='A':
                    match_gsp2.append('|')
                elif c=='R' and d =='G':
                    match_gsp2.append('|')
                elif c=='Y' and d =='T':
                    match_gsp2.append('|')
                elif c=='Y' and d =='C':
                    match_gsp2.append('|')
                elif c=='M' and d =='A':
                    match_gsp2.append('|')
                elif c=='M' and d =='C':
                    match_gsp2.append('|')

                else:
                    match_gsp2.append(' ')

            #now get the best alignments only
            all_alignments.append(alignments)
            all_matches.append(match)


            # print best_alignments[0][0]
        if best_score>64:
            outfile.write(line)
            try:
                sys.argv[2]#only print the alignments to stdout if using the -v flag, for -verbose
                print
                print "HIV"+strand,best_alignments[0]
                print '    ',"".join(match)
                print "READ",best_alignments[1]
                if strand =='+':
                    print '   '," "*(best_start+1),"".join(match_gsp2)
                    print "GSP2"," "*(best_start),gsp2_alignments[0][0]
                if strand =='-':
                    print '    ',"".join(match_gsp2)
                    print "GSP2",gsp2_alignments[0][0]
                print best_score
                print
            except IndexError:
                pass
            
outfile.close()

