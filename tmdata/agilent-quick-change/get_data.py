import os
from io import StringIO
import time
import random 
import argparse
import pandas as pd
from selenium import webdriver
from selenium.webdriver.support.ui import Select
from tqdm import tqdm

from bm3 import seq 

SEQ_OUT = 'sequences.csv'
INFO_OUT = 'info.csv'

choices = [ "A(ala)",
            "C(cys)",
            "D(asp)",
            "E(glu)",
            "F(phe)",
            "G(gly)",
            "H(his)",
            "I(ile)",
            "K(lys)",
            "L(leu)",
            "M(met)",
            "N(asn)",
            "P(pro)",
            "Q(gln)",
            "R(arg)",
            "S(ser)",
            "T(thr)",
            "V(val)",
            "W(trp)",
            "Y(tyr)"]

mutation_dropdowns_xpaths = [f'/html/body/div[1]/div/div[1]/div[2]/div/form/div[13]/div[1]/div[3]/div[{i}]/select' for i in range(2,9)] # not working
mutation_dropdowns_ids = [f'site{i}' for i in range(8,15)]
mutation_checkboxes  = [f'/html/body/div[1]/div/div[1]/div[2]/div/form/div[13]/div[1]/div[5]/div[{i}]/div' for i in range(1,494)] # 1st
design_primers_button = '//*[@id="qcpd"]/form/div[13]/div[3]/button'
clear_input_button = '/html/body/div[1]/div/div[1]/div[2]/div/form/div[15]/button'

primer_seq_table = '/html/body/div[1]/div/div[1]/div[2]/div/form/div[16]/div[2]'
oligo_info_table = '/html/body/div[1]/div/div[1]/div[2]/div/form/div[16]/div[5]'

def main(args):
    # file output - going to parse properly later
    pd.DataFrame([], columns = ["Primer Name Primer Sequence (5' to 3')"]).to_csv(SEQ_OUT)
    pd.DataFrame([], columns = ["Primer Name Length (nt.) Tm Duplex Energy at 68 °C Energy Cost of Mismatches"]).to_csv(INFO_OUT)

    # browser setup
    driver = webdriver.Chrome()
    driver.get('https://www.agilent.com/store/primerDesignProgram.jsp')

    for i in tqdm(range(args.num)):

        dna_input = driver.find_element_by_xpath('/html/body/div[1]/div/div[1]/div[2]/div/form/div[9]/div[2]/textarea')
        dna_input.send_keys(seq)
        upload_translated_button = driver.find_element_by_xpath('/html/body/div[1]/div/div[1]/div[2]/div/form/div[11]/button[2]/span/span')
        upload_translated_button.click()
        select_mutate_checkbox = driver.find_element_by_xpath('/html/body/div[1]/div/div[1]/div[2]/div/form/div[13]/div[1]/div[2]/strong/label/input[1]')
        select_mutate_checkbox.click()

        for j in mutation_dropdowns_ids:
            select = Select(driver.find_element_by_id(j))
            select.select_by_visible_text(random.choice(choices))
            # mutations not allowed near ends
            check = driver.find_element_by_xpath(random.choice(mutation_checkboxes[5:-5])) 
            check.click()
            
        button = driver.find_element_by_xpath(design_primers_button)
        button.click()

        seq_table = driver.find_element_by_xpath(primer_seq_table)
        info_table = driver.find_element_by_xpath(oligo_info_table)
        df_seq = pd.read_csv(StringIO(seq_table.text))
        df_info = pd.read_csv(StringIO(info_table.text))
        df_seq.to_csv(SEQ_OUT, mode = 'a', header = False, index = False)
        df_info.to_csv(INFO_OUT, mode = 'a', header = False, index = False)

        button = driver.find_element_by_xpath(clear_input_button)
        button.click() # takes longer, saves me unchecking all checkboxes

        



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--num', type = int, default = 5)
    args = parser.parse_args()
    main(args)
