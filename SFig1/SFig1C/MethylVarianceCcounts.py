import numpy as np
import os,sys

def main():
	input_info = open(sys.argv[1],"r")
	output_info = open(sys.argv[2],"w")
	id_column = int(sys.argv[3])-1

	methyl_dic = {}

	for each in input_info:
		each = each.strip().split()
		region_id = each[id_column]
		methylratio = float(each[3])/float(each[4])

		if region_id in methyl_dic:
			methyl_dic[region_id].append(methylratio)
		else:
			methyl_dic[region_id] = [methylratio]


	for each in methyl_dic:
		c_counts = len(methyl_dic[each])
		var = np.var(methyl_dic[each])
		std = np.std(methyl_dic[each])
		mean = np.mean(methyl_dic[each])
		output_info.write("{}\t{}\t{}\t{}\t{}\n".format(each,var,std,mean,c_counts))

	input_info.close()
	output_info.close()

main()