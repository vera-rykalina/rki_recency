# Import libraries
import pandas as pd
import sys
import re

infilename = sys.argv[1]
outfilename = sys.argv[2]


# Read .csv file
f = open(infilename, "r")
df = pd.read_csv(f, sep = ",")
f.close()


name1 = infilename.rsplit("/")[-1] # gives a file name.csv
name2 = name1.split("_BaseFreqs_")[0] # gives a sample ID
#name3 = name1.split("_Rega_")[-1].split(".")[-2] # {run_index}_{framgment}_20M


# Select only what is needed
#df = df.loc[:,["name", "assignment", "pure", "crf"]]

# Find Max (throughout columns A count: T count)
df["Max"] = df[["A count", "C count", "G count", "T count"]].max(axis=1)

# Find Sum (throughout columns A count: T count)
df["Sum"] = df[["A count", "C count", "G count", "T count"]].sum(axis=1)

# Find Max (throughout columns A count: T count)
df["MAF"] = 1 - df["Max"]/df["Sum"]


# Prepare a .csv file
df.to_csv(name2 + ".csv", sep=",", index=False, encoding="utf-8")

print(df.head())