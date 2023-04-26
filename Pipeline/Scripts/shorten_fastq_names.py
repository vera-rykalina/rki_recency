# Import libraries
import sys
import re


infilename = sys.argv[1]
outfilename = sys.argv[2]


pattern = r"\d{6}_\d{2}-\d{5}_(HIV\d{2}-\d{5})_\d{2}_\d{3}_\w{3}_L[0-9]{3}_(R\d)_\d{3}.(fastq.gz)"
substring = r"\1_\2.\3"

#names = ""

for infilename in sys.argv[1:]:
    outfilename = re.sub(pattern, substring, infilename)



