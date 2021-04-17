from io import StringIO
import re
import time
import pandas as pd
import swalign as sw
from tqdm import tqdm
from bm3 import seq as SEQ

def clean_data(seq_path, info_path):
    sequences = pd.read_csv(seq_path, header=None).iloc[:,0].str.extract('\-([acgt]+)\-')
    info = pd.read_csv(info_path, header=None).iloc[:,0].str.extract('(\d+\.\d+)℃')
    df = pd.concat([sequences, info], axis=1).dropna()
    df.columns = ['seq','tm']
    return df


def rev_comp(s):
    comp = {'a':'t','t':'a','c':'g','g':'c'}
    return ''.join([comp[i] for i in s[::-1]])


REV = rev_comp(SEQ)

def split_primer(s):
    # split into left flank, mismatch, right flank
    scoring = sw.NucleotideScoringMatrix(match = 2 , mismatch = -1)

    aligner = sw.LocalAlignment(scoring)
    aln1 = aligner.align(s, SEQ)
    aln2 = aligner.align(s, REV)
    if aln1.identity > aln2.identity:
        best_aln = aln1 
    else:
        best_aln = aln2
    f = StringIO()
    best_aln.dump(out = f)
    lines = f.getvalue().split('\n')
    # 0th line is query, 2nd line is ref 
    query_aln = re.findall('\d ([actg-]+) \d', lines[0])[0]
    ref_aln = re.findall('\d ([actg-]+) \d', lines[2])[0]
    # find diff and split 
    mismatches = [i for i, (j,k) in enumerate(zip(ref_aln, query_aln)) if j != k]
    left_flank = query_aln[:min(mismatches)]
    if len(mismatches) > 1:
        mismatch = query_aln[min(mismatches):max(mismatches)]
    else:
        mismatch = query_aln[mismatches[0]]
    right_flank = query_aln[max(mismatches):]
    return left_flank.replace('-',''), mismatch.replace('-',''), right_flank.replace('-','')




def main():
    df = clean_data('sequences.csv', 'info.csv')

    with open('clean_data.csv','a') as f:
        f.write('l,m,r,tm\n')
        for i, t in tqdm(zip(df['seq'], df['tm']), total = len(df)):
            l, m, r = split_primer(i)
            f.write(f'{l},{m},{r},{t}\n')

if __name__ == '__main__':
    main()
