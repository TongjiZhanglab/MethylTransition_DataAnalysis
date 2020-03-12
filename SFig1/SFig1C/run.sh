# awk -v "OFS=\t" '{printf("%s\t%d\t%d\t%s\t%s\t%s\n",$1,($2+$3)/2-2000,($2+$3)/2+2000,$4,$5,$6)}' ~/Data/hg19/hg19.refseq.bed | sort -k1,1 -k2,2n > hg19.refseq.GenebodyRegion.sorted
# awk -v "OFS=\t" '{printf("%s\t%d\t%d\t%s\t%s\t%s\n",$1,($2+$3)/2-2000,($2+$3)/2+2000,$4,$5,$6)}' ~/Data/hg19/hg19_CGI_c6.bed | sort -k1,1 -k2,2n | grep -v - > hg19.refseq.CGIRegion.sorted
# PromoterRegion=/mnt/Storage/home/zhaochengchen/Data/hg19/hg19.refseq.Promoter.sorted
# GenebodyRegion=hg19.refseq.GenebodyRegion.sorted
# CGIRegion=hg19.refseq.CGIRegion.sorted
# CGI=/mnt/Storage/home/zhaochengchen/Data/hg19/hg19_CGI_c6.bed
# RandomRegion=/mnt/Storage/home/zhaochengchen/Data/hg19/hg19.RandomRegion.bed
# DNAmethyl=scBS_2C_11_2.bed

for sample in scBS_2C_6_1 scBS_2C_6_2 scBS_4C_7_1 scBS_4C_7_2 scBS_4C_7_3 scBS_4C_7_4 scBS_8C_9_1 scBS_8C_9_2 scBS_8C_9_3 scBS_8C_9_4 scBS_8C_9_5 scBS_8C_9_6 scBS_8C_9_7 scBS_8C_9_8
do
	# for region in PromoterRegion GenebodyRegion CGIRegion CGI RandomRegion
	# do
	# 	intersectBed -wo -a ${sample}.bed -b ${!region} > DNAmethyl_${sample}_${region}.txt
	# done

	# python MethylVarianceCcounts.py DNAmethyl_${sample}_PromoterRegion.txt DNAmethylVar_${sample}_PromoterRegion.txt 9
	# python MethylVarianceCcounts.py DNAmethyl_${sample}_GenebodyRegion.txt DNAmethylVar_${sample}_GenebodyRegion.txt 9
	# python MethylVarianceCcounts.py DNAmethyl_${sample}_CGIRegion.txt DNAmethylVar_${sample}_CGIRegion.txt 10
	# python MethylVarianceCcounts.py DNAmethyl_${sample}_CGI.txt DNAmethylVar_${sample}_CGI.txt 10
	# python MethylVarianceCcounts.py DNAmethyl_${sample}_RandomRegion.txt DNAmethylVar_${sample}_RandomRegion.txt 9
	Rscript SFig1C.DNAmethylVar.r ${sample}
done