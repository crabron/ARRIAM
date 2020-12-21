#! /usr/bin/env python

import pandas
import numpy as np
import os

df = pandas.read_table('~/storage/mal/texmap.txt', sep='\t', index_col=0)
os.makedirs("script")
for tax, values in df.iterrows():
    row1 = df.loc[tax]
    df_i = pandas.DataFrame()
    for name1, values1 in row1.iteritems():
        l1 = []
        l2 = []
        for name2, values2 in row1.iteritems():
            if values2 == 0:
                a = "NaN"
            else:    
                a = values1/values2
            l1.append(name2)
            l2.append(a)
        s1 = pandas.Series(l2, index=l1)
        df_i = pandas.concat([df_i, s1], axis=1)
    df_i.columns = l1
    with open("script/{}.tsv".format(tax), "w") as oppai:
        df_i.to_csv( oppai , sep='\t')

