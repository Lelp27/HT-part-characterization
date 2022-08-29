"""
1. Seq 읽기
2. [:200, :-200]

양끝 200씩만 가지고 Adapter 분석 + Porechop 금지

"""
from multiprocessing import Process, Manager
from time import time
from pathlib import Path
import argparse
import gzip
from Bio import SeqIO, pairwise2
import pandas as pd
import logging
import mappy
import sys
import numpy as np
from io import StringIO


# Sample 210701 DH5a_Plib

def get_args():
    parser = argparse.ArgumentParser()
    #그룹은 둘중에 하나 필수 지정 안된값은 None
    parser.add_argument('-i', type=str, help="input Path of fastq folder")
    parser.add_argument('--tag', type=str, help="Tags fasta file", default = None)
    parser.add_argument('-o', type=str, default = None, help="")
    parser.add_argument('--format', type=str, choices = ['fastq.gz','fastq'], help="file types", default = "fastq.gz")
    parser.add_argument('-p', type=float, help="Remove parameter")
    parser.add_argument('-t', type=int, default = 8, help="multi-processing thread")
    args = parser.parse_args()

    return (args)

def load_tag_from_df(df):
    fasta = '\n'.join(['>' + str(name) + '\n' + str(seq) for name, seq in zip(df.index, df['Seq'])])
    return ([i for i in SeqIO.parse(StringIO(fasta), 'fasta')])

def load_fastq(input):
    seq_list = [record for record in SeqIO.parse(input,'fastq')]
    return (seq_list)

def multi_align(query_list, DB, output):
    for query in query_list:
        alignment = pairwise2.align.localms(query.seq, DB, 3, -6, -5, -2, one_alignment_only = True)[0]
        output.append([alignment.score, alignment.start, alignment.end])

def Align(query, DB):
    return (pairwise2.align.localms(query, DB, 3, -6, -5, -2, score_only = True))

# Input = Specific fastq file. (One or Multiple)
input_path = '/mnt/ont/2021/210701_DH5a-PLib/dh5a-plib/20210701_0930_MC-110372_0_agw324_ac3da2d2/fastq_pass/barcode01/single_read/filter.fastq.gz'
tag_db = './data/tag_db.xlsx'
args = argparse.Namespace(i = input_path, tag = tag_db, format = 'fastq.gz', t=8)
manager = Manager()

# Logging information
logging.basicConfig(
    #filename='example.log',
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S'
)
logging.INFO('RUN tag demultiplexing')
logging.INFO(f'Tag DB : {args.tag}')

# Load Sequence
if args.format == 'fastq.gz':
    seq_list = load_fastq(gzip.open(args.i, 'rt'))
elif args.format == 'fastq':
    seq_list = load_fastq(args.i)
logging.INFO(f'Loading {args.format} data ..')
logging.INFO(f'{len(seq_list)} are loadded.')

# Make Tag directory
Path.mkdir(Path(args.i).parent / Path('tag'), mode=777, exist_ok=True)
logging.INFO(f'Mkdir demultiplexed read directory : {Path(args.i).parent / Path("tag")}')

# Load tag DB
forward_tag = load_tag_from_df(pd.read_excel(args.tag, sheet_name='Forward', index_col=0))
reverse_tag = load_tag_from_df(pd.read_excel(args.tag, sheet_name='Reverse', index_col=0))

def Adapter_strand_scan(query, DB, threshold): #query = seq.record
    outlist = [[] for i in range(len(DB)+1)]
    for i in range(len(query)):
        score = np.zeros(len(DB))
        for j in range(len(DB)):
            if Align(query[i].seq, DB[j]) > len(DB[j])*threshold:
                score[j] = Align(query[i].seq, DB[j])
        if score.max() == 0:
            outlist[-1].append(i)
        else:
            outlist[score.argmax()].append(i)
    return outlist

result = Adapter_strand_scan(seq_list, DB = [forward_tag[0].seq, reverse_tag[0].seq], threshold = 1.4)

Adapter = [forward_tag[0].seq, reverse_tag[0].seq.reverse_complement()]
def dual_adapter_scan(seq_list, Adapter = []):
    pass

outlist = []
for query in seq_list:
    tmp_list = []
    for adapter in Adapter:
        plus = Align(query.seq[:200], adapter)
        minus = Align(query.seq[-200:], adapter)
        if max(plus, minus) < len(adapter)*3*0.6:
            outlist.append(0)
            break
        if plus > minus:
            tmp_list.append(1)
        elif plus < minus:
            tmp_list.append(-1)

    #print (tmp_list)
    if tmp_list == []:
        continue
    if sum(tmp_list) == 0:
        outlist.append(tmp_list[0])

n = 3
Align(seq_list[n].seq, forward_tag[0].seq)
Align(seq_list[n].seq, forward_tag[0].seq.reverse_complement())
Align(seq_list[n].seq, reverse_tag[0].seq)
Align(seq_list[n].seq, reverse_tag[0].seq.reverse_complement())

tmp = [reverse_tag[0].seq.reverse_complement() for _ in range(1000)]

Align(seq_list[0].seq[:200], forward_tag[0])
Align(seq_list[0].seq[-200:], forward_tag[0])

seq_list[0].seq[:200]


len(seq_list[0].seq[-200:])

for i in range(0, len(seq_list), 1000):
    A = Align

# 46.8
# threshold = 
start = time()

[Align(i.seq[:200], forward_tag[0]) for i in seq_list]
[Align(i.seq[:-200], forward_tag[0]) for i in seq_list]

A = [Align(i.seq[:200], forward_tag[0]) for i in seq_list[:1000]]
print (time() - start)

Align(seq_list[3].seq[:200], forward_tag[0])


final_list = manager.list()
process1 = Process(target=multi_align, args=[seq_list[:1000], forward_tag[0].seq, final_list])
process1.start()
process1.join()

print (final_list)

n = 0
for i in final_list:
    if i[0] > 45:
        n += 1






# = 전부 ALign 하기엔 속도가 너무 부족
#  Score only로 따버린 후
#  Start position 분석을 추후에 해야할 듯 ?


def Adapter_strand_scan(query, DB, threshold): #query = seq.record
    outlist = [[] for i in range(len(DB)+1)]
    for i in range(len(query)):
        score = np.zeros(len(DB))
        for j in range(len(DB)):
            if Align(query[i].seq, DB[j]) > len(DB[j])*threshold:
                score[j] = Align(query[i].seq, DB[j])
        if score.max() == 0:
            outlist[-1].append(i)
        else:
            outlist[score.argmax()].append(i)
    return outlist

result = Adapter_strand_scan(seq_list, [for_adapter], 1.4)