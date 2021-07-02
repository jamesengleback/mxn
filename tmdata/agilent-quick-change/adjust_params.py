import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import primer3

def scale(params):
    mv_conc, dv_conc, dntp_conc, dna_conc, dmsoPerc = params
    # scale - 
    # mv_conc 0 - 100 mM (monovalent cation)
    # dv_conc 0 - 100 mM (divalent cation)
    # dntp_conc 0 - 1? - check kit - mM
    # dna_conc ~ 50 nm
    #mv_conc *= 100
    #dv_conc *= 100
    #dntp_conc *= 1
    #dna_conc *= 100
    return mv_conc, dv_conc, dntp_conc, dna_conc, dmsoPerc


def calcTm(l,m,r, params):
    # params is torch tensor
    mv_conc, dv_conc, dntp_conc, dna_conc, dmsoPerc = scale(params)
    ltm = primer3.bindings.calcTm(l, 
            mv_conc = mv_conc, 
            dv_conc = dv_conc, 
            dntp_conc = dntp_conc, 
            dna_conc = dna_conc,
            tm_method = 'santalucia',
            salt_corrections_method = 'santalucia')
    rtm = primer3.bindings.calcTm(r, 
            mv_conc = mv_conc, 
            dv_conc = dv_conc, 
            dntp_conc = dntp_conc, 
            dna_conc = dna_conc,
            tm_method = 'santalucia',
            salt_corrections_method = 'santalucia')
    m_hairpin = primer3.bindings.calcHairpin(m, 
            mv_conc = mv_conc, 
            dv_conc = dv_conc, 
            dntp_conc = dntp_conc, 
            dna_conc = dna_conc)
    mtm = m_hairpin.tm
    return sum([ltm, rtm]) - 0.75 * dmsoPerc

def main():
    df = pd.read_csv('clean_data.csv')

    # these work ok
    params = [25, # mv - 25-100
            0.5, # dv - 0.5 - 5
            0.8, # dntp
            1, # dna
            8] # dmso < 10%
    y = []
    yh = []
    for l,m,r, tm in zip(df['l'], df['m'], df['r'], df['tm']):
        tm_pred = calcTm(l,m,r,params)
        y.append(tm)
        yh.append(tm_pred)

    sns.kdeplot(y,yh, fill = True)
    plt.scatter(y,yh, marker='+', alpha=0.5)
    plt.plot([0,100],[0,100])
    plt.xlim(0,100)
    #plt.ylim(0, 120)
    plt.xlabel('true')
    plt.ylabel('pred')
    plt.savefig('params-accuracy.png')
    plt.show()

if __name__ == '__main__':
    main()
