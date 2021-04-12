import argparse
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split
from tqdm import tqdm
import matplotlib.pyplot as plt

class Data(Dataset):
    def __init__(self, path, nrows = None):
        self.df = pd.read_csv(path, nrows = nrows)
        self.max_seq_len = max(self.df['seq'].str.len())
        self.x = torch.cat([self.ohe(i).unsqueeze(0) for i in self.df['seq']], dim = 0)
        self.y = torch.from_numpy(self.df[['tm']].values).float() 
        self.y -= self.y.min()
        self.y /= self.y.max()
    def __getitem__(self, idx):
        return self.x[idx], self.y[idx]
    def __len__(self):
        return len(self.df)
    def ohe(self, seq):
        x = torch.zeros((4, self.max_seq_len))
        d = dict(zip('ACGT',range(4)))
        for i, j in enumerate(seq):
            x[d[j],i] = 1.
        return x
            

class Model(nn.Module):
    def __init__(self):
        super().__init__()
        self.c1 = nn.Conv1d(in_channels = 4, out_channels = 64, kernel_size = 3, padding = 1, bias = False)
        self.c2 = nn.Conv1d(in_channels = 64, out_channels = 32, kernel_size = 5, padding = 2, bias = False)
        self.c3 = nn.Conv1d(in_channels = 32, out_channels = 16, kernel_size = 5, padding = 2, bias = False)
        self.c4 = nn.Conv1d(in_channels = 16, out_channels = 8, kernel_size = 5, padding = 2, bias = False)
        self.l1 = nn.Linear(in_features = 240, out_features = 1, bias = False)

    def forward(self, x):
        batch_size = x.shape[0]
        x = torch.relu(self.c1(x))
        x = torch.relu(self.c2(x))
        x = torch.relu(self.c3(x))
        x = torch.relu(self.c4(x))
        x = x.reshape(batch_size, -1)
        return torch.relu(self.l1(x))

# class Model(nn.Module):
#     def __init__(self):
#         super().__init__()
#         self.c1 = nn.Conv1d(in_channels = 4, out_channels = 64, kernel_size = 3, padding = 1)
#         self.c2 = nn.Conv1d(in_channels = 64, out_channels = 32, kernel_size = 5, padding = 2)
#         self.lstm = nn.LSTM(input_size = 32, hidden_size = 32, bidirectional=True)
#         self.l1 = nn.Linear(in_features = 64, out_features = 1)
# 
#     def forward(self, x):
#         batch_size = x.shape[0]
#         x = self.c1(x)
#         x = torch.relu(x)
#         x = self.c2(x)
#         x = torch.relu(x)
#         x = x.permute(2,0,1)
#         out, (h0, c0) = self.lstm(x)
#         h0 = h0.permute(1,2,0)
#         h0 = h0.reshape(batch_size, -1)
#         h0 = torch.relu(h0)
#         yh = self.l1(h0)
#         yh = torch.sigmoid(yh)
#         return yh

def main(args):
    data = Data('q5tm.csv', nrows = args.nrows)
    if args.cuda:
        data.cuda()
    trainSize = round(0.8*len(data))
    testSize = len(data) - trainSize
    train, test = random_split(data, [trainSize, testSize])
    trainLoader = DataLoader(train, batch_size = 32, shuffle = True, num_workers = 2)
    testLoader = DataLoader(test, batch_size = 64, shuffle = True, num_workers = 2)

    model = Model()
    opt = torch.optim.Adam(model.parameters(), lr = 1e-4)
    loss_fn = nn.MSELoss()

    if args.cuda:
        model.cuda()

    if args.wandb:
        wandb.init('mxn-q5')

    for epoch in range(10):
        for xi, yi in tqdm(trainLoader):
            yh = model(xi)
            loss = loss_fn(yh, yi)
            loss.backward()
            opt.step()
            opt.zero_grad()
            wandb.log({'loss':loss})
        print(loss)

    model.eval()
    yh = torch.tensor([])
    y = torch.tensor([])

    for xi, yi in tqdm(testLoader):
        yhi = model(xi)
        yh = torch.cat((yh, yhi.detach()))
        y = torch.cat((y, yi.detach()))

    # descale yh
    #yh *= data.df['tm'].max()
    #yh += data.df['tm'].min() 

    ## descale y
    #y *= data.df['tm'].max()
    #y += data.df['tm'].min() 

    plt.scatter(yh.numpy(),y.numpy(), s = 0.4)
    plt.plot([0,1],[0,1], lw = 0.2)
    plt.xlabel('pred °C')
    plt.ylabel('actual °C')
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--nrows', type=int)
    parser.add_argument('-w','--wandb', action='store_true')
    parser.add_argument('-c','--cuda', action='store_true')
    args = parser.parse_args()
    main(args)
