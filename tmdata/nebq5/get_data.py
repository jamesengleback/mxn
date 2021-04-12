import re
import os
import random
import time
import argparse
from tqdm import tqdm
import pandas as pd
from selenium import webdriver

def random_primer(minLen, maxLen):
    l = random.randint(minLen, maxLen)
    return ''.join(random.choices('ACTG',k=l))

def main(args):
    driver = webdriver.Chrome()
    driver.get('http://tmcalculator.neb.com/#!/main')
    
    time.sleep(0.1)
    
    input_1 = driver.find_element_by_xpath('//*[@id="p1"]')
    input_2 = driver.find_element_by_xpath('//*[@id="p2"]')
    p1_tm = driver.find_element_by_xpath('//*[@id="tm1"]/div[2]/strong[3]')
    p2_tm = driver.find_element_by_xpath('//*[@id="tm2"]/div[2]/strong[3]')
    anneal_at = driver.find_element_by_xpath('//*[@id="ta"]/h2')


    save_path = 'q5tm.csv'
    if not os.path.exists(save_path):
        pd.DataFrame([], columns = ['seq','tm']).to_csv(save_path,index=False)

    for i in tqdm(range(args.num)):
        p1 = random_primer(args.lenMin, args.lenMax)
        p2 = random_primer(args.lenMin, args.lenMax)
        input_1.send_keys(p1)
        input_2.send_keys(p2)


        try:
            t1 = float(re.findall('\d+', p1_tm.text)[0])
            t2 = float(re.findall('\d+', p1_tm.text)[0])
            #t_comb = float(re.findall('\d+', anneal_at.text)[0])
            df = pd.DataFrame([[p1, t1],[p2,t2]], columns = ['seq','tm'])
            df.to_csv(save_path, mode='a',header=None, index = False)
        except:
            pass

        input_1.clear()
        input_2.clear()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--lenMin',type=int, default=8)
    parser.add_argument('-j', '--lenMax',type=int, default = 30)
    parser.add_argument('-n', '--num',type=int, default = 5)
    args = parser.parse_args()
    main(args)

