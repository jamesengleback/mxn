import mxn
from bm3 import bm3, orf

def main():
    with open('mutations','r') as f:
        mutations = f.read().splitlines()

    cds = mxn.CDS(orf, bm3)
    cds.mutate(88,'A')
    cds.mutate(33,'P')
    print(cds.mutations)

if __name__ == '__main__':
    main()
