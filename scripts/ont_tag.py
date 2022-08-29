import argparse, os, sys
from glob import glob
from Bio import SeqIO, pairwise2
import numpy as np
from time import time


def main():
    parser = argparse.ArgumentParser()
    #그룹은 둘중에 하나 필수 지정 안된값은 None
    parser.add_argument('-i', type=str, help="input Path of fastq folder")
    parser.add_argument('-b', type=str, help="barcode fasta file", default = None)
    parser.add_argument('-o', type=str, default = None, help="")
    parser.add_argument('-t', type=str, choices = ['fasta','fastq'], help="file types", default = "fastq")
    
    parser.add_argument('-s', action='store_true', help="Save background removed image files")

    global args
    args = parser.parse_args()

def rev_com(seq):  #rev_com = reverse complementary함수.
    seq = seq[::-1]
    table = str.maketrans('ATGC','TACG')
    return seq.translate(table)

def align(query, DB):
    #localds = PWM 적용가능
    alignments = pairwise2.align.localms(query, DB, 2, -2, -3, -2, score_only = True)
    return alignments

def Adapter_strand_scan(query, DB, threshold): #query = seq.record
    outlist = [[] for i in range(len(DB)+1)]
    for i in range(len(query)):
        score = np.zeros(len(DB))
        for j in range(len(DB)):
            if align(query[i].seq, DB[j]) > len(DB[j])*threshold:
                score[j] = align(query[i].seq, DB[j])
        if score.max() == 0:
            outlist[-1].append(i)
        else:
            outlist[score.argmax()].append(i)
    return outlist

def tag_scan(query, tag_DB, align_threshold=1.4, tag_threshold=0.1, PWM = None):
    tag_list = [[] for i in range(len(tag_DB))]
    for i in range(len(query)):
        score = np.zeros(len(tag_DB))
        for j in range(len(tag_DB)):
            if align(query[i].seq, tag_DB[j].seq) > len(tag_DB[j])*align_threshold:
                score[j] = align(query[i].seq, tag_DB[j].seq)/(len(tag_DB[j])*2)
        tag_score = score.max() - np.partition(score, -2)[-2]
        if tag_score > tag_threshold:
            tag_list[score.argmax()].append(i)
    return tag_list

if __name__ == '__main__':
    main()

#Primer는 Fasta file을 읽자.
if args.o == None:
    args.o = args.i
primer = r'C:\Users\user\task\Row_data\DB\M13_tag.fasta'
M13F = 'GTAAAACGACGGCCAGT'

seq_list = glob(args.i + '\*.fastq')

#tag폴더생성
try:
    if not(os.path.isdir(args.i + '\\tag')):
        os.makedirs(args.i +'\\tag')
except:
    print ("Tag Directory already is")
    if input("Overwrite Tag ? [Y/N]") == "Y":
        pass
    else:
        sys.exit("Break! tag folder is already in")

#----------
summary = open(args.i + '\\summary.txt', 'w')
summary.writelines('Barcode Path = ' + args.i + '\n')
summary.writelines('Primer Path = ' + primer + '\n')
#----------

#모든 sequence parsing
library = []
for i in seq_list:
    [library.append(record) for record in SeqIO.parse(i, args.t)]

tag_DB = []
[tag_DB.append(record) for record in SeqIO.parse(primer, 'fasta')]
tag_DB[-1].seq
#------------
summary.writelines('전체 Read 수 : ' + str(len(library)) + '\n')
#------------

#Adapter alignment
#Adapter의 서열을 더늘려서(앞뒤로) 해상력을 높이는 것이 좋아보임.
result = Adapter_strand_scan(library, DB = [M13F, rev_com(M13F)], threshold = 1.4)

#(-)starnd to (+)
for i in result[1][::-1]:
    library[i].seq = rev_com(str(library[i].seq))
for i in result[-1][::-1]:
    library.pop(i)
#------------
summary.writelines('Adapter 발견 Read 수 : ' + str(len(library)) +'\n')
#------------

start = time()
tag_result = tag_scan(library, tag_DB, align_threshold = 1.4, tag_threshold = 0.1)
run_time = time() - start

#---------------
summary.writelines('\n' + '=======================================' + '\n')
summary.writelines("Tag_scan run time : " + str(run_time) + '\n')
for i in range(len(tag_result)):
    summary.writelines(
        tag_DB[i].id + ' : ' + str(len(tag_result[i])) + '\n'
    )
summary.writelines('unclassified : ' + str(len(library) - sum([len(i) for i in tag_result])) + '\n')
summary.writelines('\n' + '=======================================' + '\n')
summary.close()
#---------------
#save_tag
for i in range(len(tag_result)):
    ID = tag_DB[i].id
    f = open(args.i + '\\tag\\' + ID + '.fasta', 'w')
    for j in range(len(tag_result[i])):
        SEQ = library[tag_result[i][j]].seq
        f.writelines('>' + str(ID) + '\n' + str(SEQ) + '\n')
    f.close()