import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import primer3
import torch
import wandb
from tqdm import tqdm

def scale(params):
    params = params.detach().numpy()
    mv_conc, dv_conc, dntp_conc, dna_conc = params
    # scale - 
    # mv_conc 0 - 100 mM (monovalent cation)
    # dv_conc 0 - 100 mM (divalent cation)
    # dntp_conc 0 - 1? - check kit - mM
    # dna_conc ~ 50 nm
    mv_conc *= 100
    dv_conc *= 100
    dntp_conc *= 1
    dna_conc *= 100
    return mv_conc, dv_conc, dntp_conc, dna_conc  

def calcTm(l,m,r, params):
    # params is torch tensor
    mv_conc, dv_conc, dntp_conc, dna_conc = scale(params)
    ltm = primer3.bindings.calcTm(l, 
            mv_conc = mv_conc, 
            dv_conc = dv_conc, 
            dntp_conc = dntp_conc, 
            dna_conc = dna_conc)
    rtm = primer3.bindings.calcTm(r, 
            mv_conc = mv_conc, 
            dv_conc = dv_conc, 
            dntp_conc = dntp_conc, 
            dna_conc = dna_conc)
    m_hairpin = primer3.bindings.calcHairpin(m, 
            mv_conc = mv_conc, 
            dv_conc = dv_conc, 
            dntp_conc = dntp_conc, 
            dna_conc = dna_conc)
    mtm = m_hairpin.tm
    return sum([ltm, rtm])

def main():
    df = pd.read_csv('clean_data.csv')

    params = torch.randn(4, requires_grad = True)
    opt = torch.optim.Adam([params], lr = 1e-8)
    loss_fn = torch.nn.MSELoss()
    wandb.init('mxn-agilent')

    norm = lambda x : (x - df['tm'].min()) / df['tm'].max()

    for i in tqdm(range(10000)):
        y = []
        yh = []
        for l,m,r, tm in zip(df['l'], df['m'], df['r'], df['tm']):
            tm_pred = calcTm(l,m,r,torch.relu(params))
            y.append(tm)
            yh.append(tm_pred)

        loss = loss_fn(norm(torch.tensor(yh)), 
                norm(torch.tensor(y, requires_grad=True)))
        loss.backward()
        opt.step()
        opt.zero_grad()
        wandb.log({'loss':loss})

    # eval & plot
    y = []
    yh = []
    for l,m,r, tm in zip(df['l'], df['m'], df['r'], df['tm']):
        tm_pred = calcTm(l,m,r,params)
        y.append(tm)
        yh.append(tm_pred)
    plt.scatter(y,yh)
    plt.plot([0,100],[0,100])
    plt.xlim(0,100)
    plt.ylim(0, 100)
    plt.xlabel('true')
    plt.ylabel('pred')
    wandb.log({'fig':plt})

if __name__ == '__main__':
    main()
