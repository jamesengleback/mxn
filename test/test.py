import re
import pandas as pd
from tqdm import tqdm
import mxn
from bm3 import bm3, orf
pd.options.display.max_colwidth = 100

def main():
    with open('mutations','r') as f:
        mutations = f.read().splitlines()
    # parse_mutation = lambda x : re.findall('\d+', x)[0], x[-1]
    def parse_mutation(m):
        return int(re.findall('\d+', m)[0]), re.findall('\w', m)[-1] 

    cds = mxn.CDS(orf, bm3)

    for i in tqdm(mutations):
        pos, aa = parse_mutation(i)
        cds.mutate(pos, aa)
    print(pd.DataFrame(cds.primers).T)

    #print(cds.mutations)

if __name__ == '__main__':
    main()
