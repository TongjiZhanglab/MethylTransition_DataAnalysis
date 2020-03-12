from sklearn.datasets import load_boston
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import os,sys

def main():
	feature_file = sys.argv[1]
	label_file = sys.argv[2]
	out_info = open(sys.argv[3],"w")

	# X = np.transpose(np.genfromtxt("/mnt/Storage/home/zhaochengchen/Work/3.CellTracing/data/Cas9/ExpressionMatrix_e1_Layel_9_mm1.txt", delimiter="\t", skip_header=1)[:,1:])
	# Y = np.genfromtxt("/mnt/Storage/home/zhaochengchen/Work/3.CellTracing/result/RandomForest/ClassLabel_e1_Layel_9_mm1.txt", delimiter="\t", skip_header=0)[:,1]

	X = np.genfromtxt(feature_file, delimiter="\t", skip_header=1)[:,1:]
	Y = np.genfromtxt(label_file, delimiter="\t", skip_header=0)[:,1]

	rf = RandomForestRegressor()
	rf.fit(X, Y)

	for each in rf.feature_importances_:
		out_info.write(str(each)+"\n")

	out_info.close()

main()
	