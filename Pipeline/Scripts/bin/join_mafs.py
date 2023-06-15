#!/usr/local/Caskroom/miniconda/base/bin/python

# Import libraries
import pandas as pd
import sys
import functools as ft


#infilename = sys.argv[1]
#outfilename = sys.argv[2]


dfs = []
for infilename in sys.argv[1:]:
    dfs.append(pd.read_csv(infilename, sep = ","))

ref = pd.read_csv("HXB2_refdata.csv", sep = ",")
#ref = pd.DataFrame(ref["pos"], columns=["pos"])

#dfs = [df.set_index('pos') for df in dfs] 
#df = dfs[0].join(dfs[1:])
#df = ft.reduce(lambda  left, right: left.join(right, how='outer', on='pos'), dfs)
#df= ft.reduce(lambda left, right: pd.merge(left, right, on='pos'), dfs)
#df= ft.reduce(lambda left, right: pd.merge(left, right, on = "pos", how = "outer"), dfs)
#df = dfs[0]
dfs.insert(0, ref)

df=dfs[0]
print(len(dfs))
for df_ in dfs[1:]:
    df = df.merge(df_, on="pos", how="left")


df.to_csv("joined" + "_MAF" + ".csv", sep=",", index=False, encoding="utf-8")
