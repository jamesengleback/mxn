'''
primer design for mutating proteins
todo:
    primers
    score
    tmcalc.js

'''
import os
import nwalign3 as nw
from data import codons

class cds:
    def __init__(self, dna, aaseq = None):
        # dna is a single ORF
        # aaseq is the canonical aa sequence
        self.dna = dna.upper()
        self.aaseq = aaseq # reference sequence
        self.seq = self.aaseq # for changing
        self.primers = {}
    @property
    def translation(self):
        return ''.join([codons[i.upper()] for i in self.codons])
    @property
    def codons(self):
        return [self.dna[i:i+3] for i in range(0,len(self.dna),3)]
    @property
    def mutations(self):
        return utils.diff(self.aaseq, self.seq)
    def mutate(self, pos, aa):
        self.seq = utils.mxn(self.seq, pos,aa.upper())
        new_codon = utils.best_codon(codons, aa)
        c = list(self.codons)
        c[pos] = new_codon
        self.dna = ''.join(c)
        return self.seq
    def primers(self, mutation):
        # mutation codon
        # grow

        pass

class utils:
    def aln(s1,s2):
        aln1, aln2 = nw.global_align(s1,s2)
        return aln1, aln2
    def diff(s1, s2):
        return {i:{'from':x, 'to':y} for i, (x,y) in enumerate(zip(s1,s2)) if x!=y and y!='-' and x!='-'}
    def mxn(s, pos,to):
        s = list(s)
        s[pos] = to
        s = ''.join(s)
        return s
    def print_aln(aln1, aln2, name1='', name2=''):
        symbols = ''.join(['|' if i==j else ' ' for i,j in zip(aln1, aln2)])
        width = int(os.popen('stty size').read().split()[1]) - (15 + max(len(name1), len(name2)))
        for i in range(0,len(aln1), width):
            print(f"\
{i} {aln1[i:i+width]} {i+width} {name1}\n\
{' '*len(str(i))} {symbols[i:i+width]}\n\
{i} {aln2[i:i+width]} {i+width} {name2}")
    def best_codon(d,aa):
        options = {i:float(d[i]['frac']) for i in d if d[i]['aa'] == aa.upper()}
        return max(options, key=options.get)
    def score_primer(s):
        pass
def test():
    from data import bm3aa, orf
    
    c = cds(orf, bm3aa)
    s1 = c.dna
    c.mutate(100,'a')
    c.mutate(33,'e')
    aln1, aln2 = nw.global_align(orf.upper(), c.dna.upper())
    utils.print_aln(aln1, aln2,'orf','new')

if __name__ == '__main__':
    test()
