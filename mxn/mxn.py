import os
import random
import nwalign3 as nw
import primer3
import numpy as np
from scipy import optimize

from mxn.data import Codons
from mxn.agilent_tm import agilent_tm

COMPLIMENT = {'A':'T','T':'A','C':'G','G':'C'}

reverse_comp = lambda s : ''.join([d[i] for i in s[::-1]])

class CDS:
    def __init__(self, dna, aaseq = None):
        # dna is a single ORF
        # aaseq is the canonical aa sequence
        self.dna = dna.upper()
        self._dna = dna.upper() # original
        self.aaseq = aaseq # for changing
        self._aaseq = aaseq # reference sequence
        self.seq = self.aaseq # for changing - needs to be a copy
        self.primers = {}
    @property
    def codons(self):
        return [self.dna[i:i+3] for i in range(0,len(self.dna),3)]
    @property
    def translation(self):
        return ''.join([Codons[i]['aa'] for i in self.codons])
    @property
    def mutations(self):
        return utils.diff(self._aaseq, self.aaseq)
    def mutate(self, pos, aa):
        start_seq = self.dna
        self.aaseq = self.mutate_dna(pos, aa) # changes self.dna
        primers = self.make_primers(start_seq, self.dna)
        ID =f'{pos}{aa}' 
        primers['id'] = ID
        self.primers[ID] = primers
    def mutate_dna(self, pos, aa):
        self.aaseq = utils.mxn(self.aaseq, pos,aa.upper())
        new_codon = utils.best_codon(Codons, aa)
        c = list(self.codons)
        c[pos] = new_codon
        self.dna = ''.join(c)
        return self.aaseq
    def make_primers(self, seq_from, seq_to, tm=78):
        if seq_from != seq_to:
            changes_idx = [i for i, (j,k) in enumerate(zip(seq_from, seq_to)) if j != k]
            payload = seq_to[min(changes_idx):max(changes_idx)]

            left = lambda n : seq_from[max(changes_idx):max(changes_idx) + n]
            right = lambda n : seq_from[min(changes_idx) - n :min(changes_idx)]

            homoTm = lambda s : primer3.calcHomodimerTm(s,mv_conc= 25, dv_conc = 0.5)
            endscore = lambda s : sum([1 for i in s[1] + s[-2] if i == 'C' or i == 'G']\
                                    + [2 for i in s[0] + s[-1] if i == 'C' or i == 'G']) # max 6
            scorefn = lambda l, m, r : abs(agilent_tm(l, r) - tm) \
                                            - endscore(l+m+r) # + homoTm(l+m+r)
            objective = lambda ln, rn : scorefn(left(round(ln)), payload, right(round(rn)))
            helper = lambda array : objective(array[0], array[1])
            results = optimize.dual_annealing(helper, 
                                    x0 = np.array([20,20]), 
                                    bounds = ((10,40),(10,40)))
            l = left(round(results.x[0]))
            r = right(round(results.x[1]))
            score = scorefn(l,payload,r)
            tm = agilent_tm(l,r)

            primer = l + payload + r
            return {'primer':primer,
                    'tm':tm,
                    'end_score': endscore(primer),
                    'homotm':homoTm(primer),
                    'length':len(primer)}
        else:
            return {}
            
def q5tm(s):
    return primer3.calcTm(s, 
            dntp_conc = 0.2, 
            dna_conc = 1666, 
            mv_conc = 150,
            dv_conc = 0.15,
            tm_method = 'breslauer', # breslauer or santalucia
            salt_corrections_method = 'schildkraut') # schildkraut, owczarzy, santalucia

class primer:
    def __init__(self, seq):
        self.seq = seq
    def q5tm(self):
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
