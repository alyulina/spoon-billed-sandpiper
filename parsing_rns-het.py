#!/usr/bin/env python
# coding: utf-8

# parsing the output of Fedya's code for the red-necked stint to get heterozygosity averages
# averaging over the last 100 time points

import os
import numpy as np
import math

alpha = [0, -2, -4]
beta = [0, 1000, 3000, 5000, 7000]

with open('rns-het_points-to-use_fedya.txt') as f:
    points_fedya = [['mu=' + x.split()[0], 'sigma=' + x.split()[1], 'alpha=' + x.split()[2], 'beta=' + str(y)] for x in f.readlines() for y in beta]
    
rns_pop_size = 5000
n = 100 # how many time points to average over

hets_rns = []; pols_rns = []
for p in points_fedya:
        
    with open('/nfs/scistore08/kondrgrp/alyulina/sandpiper/dynamics/h=sigm/ver22/pos/rns/RNSv22_' + '_'.join(p) + '.txt') as f:
        lines_rns = f.readlines()
        
        ind = [lines_rns.index(x) for x in lines_rns if 'generation\t29999' in x] 
        lines = []
        for i in ind:
            lines.extend(lines_rns[i - n : i])
                   
        avg_het_2pq = np.mean([float(x.split()[11]) for x in lines])                                                     
        hets_rns.append(' '.join(p).replace('=', ' ') + ' het ' + '{:.12f}'.format(avg_het_2pq) + '\n')

with open('/nfs/scistore08/kondrgrp/alyulina/sandpiper/dynamics/h=sigm/ver22/pos/rns-hets-points_rns_avg-het.txt', 'w+') as o:
    o.writelines(hets_rns)

for a in alpha:
    for b in beta:
        with open('/nfs/scistore08/kondrgrp/alyulina/sandpiper/dynamics/h=sigm/ver22/pos/out/rns_points_rns_avg-het_alpha=' + '{:.1f}'.format(a) + '_beta=' + '{:.0f}'.format(b) + '.txt', 'w+') as o:
            o.writelines([x.replace('mu', '').replace('sigma', '').replace('alpha', '').replace('beta', '') for x in hets_rns if 'alpha ' + '{:.1f}'.format(a) + ' beta ' + '{:.0f}'.format(b) + ' ' in x])
