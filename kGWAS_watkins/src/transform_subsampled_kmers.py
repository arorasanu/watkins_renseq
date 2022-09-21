import sys  
import numpy as np  
import pandas as pd 

subsampled_kmer_matrix=sys.argv[1]
header_file = sys.argv[2]
output_file = sys.argv[3]

kmer_pres_list = []
with open(subsampled_kmer_matrix,'r') as f, open(header_file, 'r') as h:
	header = [l.strip() for l in h]
	for inputline in f:
		lsplit = inputline.strip().split()
		kmer_pres_list.append(list(lsplit[1]))

kmer_pres_array = np.array(kmer_pres_list).astype(int)
kmer_pres_DF = pd.DataFrame(np.transpose(kmer_pres_array), index = header)
kmer_pres_DF.to_csv(output_file, sep = "\t")