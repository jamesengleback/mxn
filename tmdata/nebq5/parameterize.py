import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import primer3
import torch
from torch.nn import MSELoss
from tqdm import tqdm

def q5tm(s, 
        dntp_conc = 0.2, 
        dna_conc = 1666, 
        dv_conc = 0.15,
        mv_conc = 0.15,
        tm_method = 'breslauer', # breslauer or santalucia
        salt_corrections_method = 'schildkraut'): # schildkraut, owczarzy, santalucia
    return primer3.calcTm(s, 
            dntp_conc = dntp_conc, 
            dna_conc = dna_conc, 
            dv_conc = dv_conc,
            tm_method = tm_method, # breslauer or santalucia
            salt_corrections_method = salt_corrections_method) # schildkraut, owczarzy, santalucia

def scaleY(y):
    return y / 100

def descaleY(y):
    return y * 100

def descaleParams(x):
    # start small torch params
    dntp_conc, dna_conc, dv_conc, mv_conc = torch.sigmoid(x).detach().numpy()
    return dntp_conc, dna_conc * 1000, dv_conc, mv_conc

def main():
    df = pd.read_csv('q5tm.csv')
    x = df['seq']
    y = df['tm']

    y_scaled = torch.tensor([scaleY(i) for i in y])

    params_out = 'params.csv'
    pd.DataFrame([], columns = ['name','dntp_conc', 'dna_conc', 'dv_conc', 'mv_conc', 'tm_method', 'salt_method', 'loss']).to_csv(params_out)

    loss_fn = MSELoss()

    batch_size = 100

    names = iter(list('ABCDEFG'))
    for tm_method in ['breslauer', 'santalucia']:
        for salt_method in ['schildkraut', 'owczarzy', 'santalucia']:

            params = torch.randn(4, requires_grad = True).float()
            opt = torch.optim.Adam([params], lr=1e-3)

            for i in tqdm(range(10000)):
                batchIdx = random.choices(x.index, k = batch_size)
                x_batch = x[batchIdx]
                y_batch = y[batchIdx]

                yh = [q5tm(i, *descaleParams(params), tm_method = tm_method, salt_corrections_method = salt_method ) for i in x_batch]
                yh_scaled = torch.tensor([scaleY(i) for i in yh], requires_grad = True)
                y_batch_scaled = torch.tensor([scaleY(i) for i in y_batch])
                loss = loss_fn(yh_scaled, y_batch_scaled)
                loss.backward()
                opt.step()
                opt.zero_grad()

            # eval and save
            yh_all = [q5tm(i, *descaleParams(params), tm_method = tm_method, salt_corrections_method = salt_method ) for i in x]
            yh_all_scaled = torch.tensor([scaleY(i) for i in yh_all], requires_grad = True)
            loss = loss_fn(yh_all_scaled, y_scaled)
            dntp_conc, dna_conc, dv_conc, mv_conc = descaleParams(params)
            name = next(names)
            pd.DataFrame([y.to_list(), yh_all], index = ['true','pred']).T.to_csv(f'{name}.csv', index=False)
            pd.DataFrame([[name, dntp_conc, dna_conc, dv_conc, mv_conc, tm_method, salt_method, loss.detach().item()]], 
                    columns = ['name','dntp_conc', 'dna_conc', 'dv_conc', 'mv_conc', 'tm_method', 'salt_method', 'loss']).to_csv(params_out, mode = 'a', header = False, index = False)
            
    fig, ax = plt.subplots(3,2, figsize = (10,10))
    for i, j in zip(ax.flatten(), 'ABCDEF'):
        data = pd.read_csv(f'{j}.csv')
        i.scatter(data['true'], data['pred'], s = 0.2)
        i.plot([0,100],[0,100])
        i.set_xlabel('pred tm °C')
        i.set_ylabel('true tm °C')
        i.set_title(j)
        i.set_ylim(0,100)
        i.set_xlim(0,100)
    plt.tight_layout()
    plt.savefig('param-fit.png')
    plt.show()


if __name__ == '__main__':
    main()
