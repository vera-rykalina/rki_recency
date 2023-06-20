#!/usr/bin/python3

# Import libraries
import pandas as pd
import sys


# Create a list of dfs
dfs = []
for infilename in sys.argv[1:]:
    dfs.append(pd.read_csv(infilename, sep = ",", index_col=False))


# Set a reference df as a df other columns should be joined
df=dfs[0]
for df_ in dfs[1:]:
    df = df.merge(df_, on="pos", how="left")

# Remove a sinity check column
df.drop(["HXB2 base"], axis=1, inplace=True)


#df = df.set_index("HXB2 base")
print(df)
# Transpose df
#df = df.T
df = df.T.reset_index()
# Reset index
#df = df.reset_index()
print(df)

# Make pos to be a header
#headers = df.iloc[0]

#df  = pd.DataFrame(df.values[1:], columns=headers)
df = df.rename_axis(None, axis = 1)
print(df)


# Create a clean csv file
df.to_csv("joined" + "_MAF" + ".csv", sep=",", index = False, header = False, encoding="utf-8")
