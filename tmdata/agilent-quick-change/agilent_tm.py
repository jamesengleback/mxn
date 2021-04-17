import primer3

def agilent_tm(l,r):
    mv_conc, dv_conc, dntp_conc, dna_conc, dmsoPerc = [25, # mv - 25-100
            0.5, # dv - 0.5 - 5
            0.8, # dntp
            1, # dna
            8] # dmso 0=10%

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
    return sum([ltm, rtm]) - 0.75 * dmsoPerc
