#! /usr/bin/env python

from sys import argv
import numpy as np
from scipy.stats import mannwhitneyu as man
import biom
import pandas
import h5py
import re
import itertools
import sys
import os
import shutil
import traceback
import time
import collections
import argparse
import re

parser = argparse.ArgumentParser(description="Convert data to frequency, calculate which otu significant change. Use Man-Whitney stats for pvalue calculation. Run-time version of code. Input - biom table in h5py format. Output - pairwise comparison data in xlsx files(or tsv) with OTU proporion and taxa, csv log file with data of OTU quantity adjustment.Sample ID have to be in format [someIDname].[replication-identificator(number)-without-dots] Examples: C.dir23.1 or someSoil.ty!")
parser.add_argument("-t", "--treshold", help="Treshold for frequency values. Default - 0.005",action="store", dest="tresh",default='0.005')
parser.add_argument("-i", "--imput", help="imput biom table",action="store", dest="inp",required=True)
parser.add_argument("-m", "--map", help="ID mapping file. Requires text file with tab delimiter between compared samples and line break (in unix format) between sample pairs. Default - everything with everything",action='store', default="all", dest="map")
parser.add_argument("-o", "--out", help="tsv - out in tsv tables. Default - export to xlsx Excel workbook.",action='store', default="xlsx", dest="export")
parser.add_argument("-c", "--chloroplast", help="yes - delete otus there identificated as chloroplasts ",action='store', default="no", dest="chloroplast")


args = parser.parse_args()


chloroplast = args.chloroplast
tresh = args.tresh
export = args.export
inp = args.inp
if args.map == 'all':
    otu_table = biom.load_table(inp)
    sample_list = otu_table.ids(axis='sample')
    atata = []
    for a in sample_list:
        i = a.split('.')
        ca = i[:-1]
        caa = ".".join(ca)
        atata.append(caa)
    at_set = set(atata)
    w=[i for i in itertools.combinations(at_set, 2)]
    with open('temp2.txt', 'a') as ID:
        for i in w:
            ID.write(i[0] + '\t' + i[1]+ '\n')
    ID = 'temp2.txt'
else:
    ID = args.map

def benchmark(func):
    def wrapper(*args, **kwargs):
        t = time.clock()
        res = func(*args, **kwargs)
        print(func.__name__, time.clock() - t)
        return res
    return wrapper

def del_singltones(inp):
    otu_table = biom.load_table(inp)
    a = [i for i in otu_table.iter_data(axis="observation")]
    b = [i for i in a if sum(i)>5]

def ratiometr(ftable, ilist):  
    left = [i for i in ftable.columns.values if i.rsplit('.', 1)[0] == ilist[0]]
    right = [i for i in ftable.columns.values if i.rsplit('.', 1)[0] == ilist[1]]
    both = left + right
    summ = ftable[both].sum()
    rtable = ftable[both]/summ
    table = pandas.concat([rtable, ftable.taxonomy], axis=1, join='inner')
    return left, right, table, summ

def ratfilt(left, right, table, tresh):
    table['mean1'] = table[left].apply(lambda row: sum(row)/len(row),axis=1)
    table['mean2'] = table[right].apply(lambda row: sum(row)/len(row),axis=1)
    ftable = table.ix[(table['mean1'] >= float(tresh)) | (table['mean2'] >= float(tresh))]
    return ftable


def man_y(rttable,left, right):
    def many(x,y):
        f = man(x,y, use_continuity=True, alternative=None)[1]
        return f
    
    rttable['pvalue'] = rttable.apply(lambda x: many(x[left],x[right]),axis=1)
    return rttable

def filter_otu(inp, idfs):
    '''
    Filter biom table by sample Id from mapping txt file.
    '''
    otu_table = biom.load_table(inp)
    slist = otu_table.ids(axis='sample')
    ll = [i for i in slist if i.startswith(idfs[0]) or i.startswith(idfs[1])]
    new_table = otu_table.filter(ll,axis='sample', inplace=False)
    return new_table

def log(out, table,ilist):
    idfs = '_'.join(ilist)
    summ = table.count()[0]
    t_unic1 = table[(table.mean1 > 0) & (table.mean2 == 0)].count()[0]
    t_unic2 = table[(table.mean2 > 0) & (table.mean1 == 0)].count()[0]
    table_gen = table.ix[(table.mean1 > 0) & (table.mean2 > 0)]
    general = table_gen.count()[0]
    increase_in_1 = table_gen[(table_gen.mean1 > table_gen.mean2) & (table_gen.pvalue<= 0.05)].count()[0]
    decrease_in_1 = table_gen[(table_gen.mean1 < table_gen.mean2) & (table_gen.pvalue<= 0.05)].count()[0]
    same = table_gen[table.pvalue > 0.05].count()[0]
    a = [summ, t_unic1, t_unic2, general ,increase_in_1, decrease_in_1, same]
    out.loc[idfs] = a
    return out, table_gen

def out(ilist, table, left, right, table_gen, export):
    a = "_".join(ilist)
    t_unic1 = table.ix[(table.mean1 > 0) & (table.mean2 == 0)]
    df2 = t_unic1[['taxonomy','mean1','mean2']]
    t_unic2 = table.ix[(table.mean2 > 0) & (table.mean1 == 0)]
    df3 = t_unic2[['taxonomy','mean1','mean2']]
    table_gen_change = table_gen.ix[table_gen.pvalue <= 0.05]
    table_gen_change['proportion'] = table_gen_change.apply(lambda x: x['mean1']/x['mean2'],axis=1)
    df1 = table_gen_change[['taxonomy','mean1','mean2','proportion']]
    if export == 'xlsx':
        writer = pandas.ExcelWriter('script/{}.xlsx'.format(a))
        df1.to_excel(writer,'general_with_proportion')
        df2.to_excel(writer,'unic_for_1')
        df3.to_excel(writer,'unic_for_2')
        writer.save()
    else:
        with open("script/{}_general.csv".format(a), "w") as sd:
            df1.to_dense().to_csv( sd , sep='\t', index = False)
        with open("script/{}unic_for_1.csv".format(a), "w") as sf:
            df2.to_dense().to_csv( sf, sep='\t', index = False)
        with open("script/{}_unic_for_2.csv".format(a), "w") as sg:
            df3.to_dense().to_csv( sg, sep='\t',  index = False)

    return

def transform(table):
    out = table.to_dataframe()
    return out

def alltransform(table):
    ftable = biom.load_table(table)
    d = []
    for i,e,w in ftable.iter(dense=True, axis='observation'):
        for a in w.values():
            # ar = a[::-1]
            vall = "|".join(a)
            d.append(vall)
    out = ftable.to_dataframe()
    out['taxonomy'] = d
    return out

def pfilter(table, ilist):
    ids = [i for i in table.columns.values if i.rsplit('.', 1)[0] == ilist[0] or i.rsplit('.', 1)[0] == ilist[1]]
    ids.append('taxonomy')
    ftable = table[ids]
    return(ftable)

def chlfilt(df, chl):
    if chl == 'yes':
        table = df[df.taxonomy.str.contains('Chloroplast') == False]
    else:
        table = df
    return table

def main(inp, ID):
    try:
        if not os.path.exists("script"):
            os.makedirs("script")
        d_out = pandas.DataFrame(columns=('summ' ,'unic_for_1',  'unic_for_2', 'general', 'increase_in_1', 'decrease_in_1', 'same'))
        with open(ID, "r") as idfile:
            idf = idfile.readlines()
            lenidf = len(idf)
            if lenidf >= 0:
                table = alltransform(inp)
            fjl=0
            ss = len(idf)
            for i in idf:

                sys.stdout.write('\r')
                fjl+=1
                sys.stdout.write("{} {}/{}".format(i,fjl,ss))
                sys.stdout.flush()

                nonn_line = i.rstrip()
                ilist = nonn_line.split("\t")
                if lenidf <= 0:
                    ftable = filter_otu(inp, ilist)
                    datatable = transform(ftable)
                else:
                    datatable = pfilter(table, ilist)
                datatable = chlfilt(datatable, chloroplast)
                left, right, rttable, summ = ratiometr(datatable, ilist)
                frttable = ratfilt(left, right, rttable, tresh)
                frttable = man_y(frttable,left, right)
                global table_out
                table_out, table_gen = log(d_out, frttable, ilist)
                out(ilist, frttable,left, right, table_gen, export)
        l = open("script/log.csv", "w")
        table_out.to_csv( l , sep='\t')
        if os.path.isfile("temp2.txt"):   
            os.remove('temp2.txt')

    except Exception:
        if os.path.isfile("temp2.txt"):   
            os.remove('temp2.txt')
        if os.path.exists("script"):
            shutil.rmtree('script', ignore_errors=True)
        print(traceback.format_exc())




if __name__ == "__main__":
    main(inp, ID)
