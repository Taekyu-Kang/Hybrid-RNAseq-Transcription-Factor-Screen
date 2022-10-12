from sys import argv
import numpy as np
import time
import pickle
import collections

#parses TF database or VCF to make nested window dictionaries

def make_TF_dict(filename):
	print "Starting to make TF dictionary..."
	start_time = time.time()
	TF = {}
	count = 0
	with open(filename) as f:
		#read in GTRD data
		for line in f:
			if line[0] == "#":
				continue
			count += 1
			linedata = line.strip().split("\t")
			chrom = linedata[0][3:]
			start = int(linedata[1])
			end = int(linedata[2])
			peak = int(linedata[3])
			name = linedata[5]

			#make new start and end locations for the binding sites +-10 bases from the ChIP-seq peak
			peak_pos = start+peak
			new_start = peak_pos-10
			new_end = peak_pos+10

			if chrom not in TF:
				TF[chrom] = []

			TF[chrom].append([new_start,new_end,name])

	print count, "TF binding sites, took", time.time()-start_time, "sec"
	for i in TF:
		print "Chromosome", i, ":", len(TF[i]), "sites."

	#start organizing the binding sites into a nested dictionary using tuples of genomic location windows as keys
	#significantly improves searching speed for future scripts

	newTF = {}

	for chrom in TF:
		if chrom in passed:
			print "Already done with chrom",chrom
			continue

		start_time = time.time()
		print "Starting to make windowed dictionary for chromosome", chrom
		newTF[chrom] = {}
		all_pos = []
		for site in TF[chrom]:
			all_pos.append(site[0])
			all_pos.append(site[1])

		maxx = max(all_pos)
		minn = min(all_pos)
		splits = (maxx-minn)/1000 #first set of windows to use as keys
		for i in range(minn,maxx,splits):
			newTF[chrom][(i,i+splits)] = {}

			splits2 = splits/10 #second set of windows to use as keys
			for j in range(i,i+splits,splits2):
				newTF[chrom][(i,i+splits)][(j,j+splits2)] = []

		newTF2 = collections.OrderedDict(sorted(newTF.items()))
		for i in newTF2:
			newTF2[i] = collections.OrderedDict(sorted(newTF2[i].items()))

		print "Done making splits, took", time.time()-start_time, "sec"

		for site in TF[chrom]:
			for window in newTF2[chrom]:
				if site[0] > window[1]:
					continue
				if site[0] in range(window[0],window[1]):
					for window2 in newTF2[chrom][window]:
						if site[0] in range(window2[0],window2[1]):
							newTF2[chrom][window][window2].append(site)
							break #breaks from the loop if nested window is found
					break

		print "Done with chromosome",chrom,"overall",time.time()-start_time,"sec"
				
		newdict = {chrom:newTF2[chrom]}
		with open("/work/tkang/primary_ASE/gtrd/bychrom/$"+chrom+"$GTRD_allTFdict_bychrom.pkl","wb") as f:
			pickle.dump(newdict,f,pickle.HIGHEST_PROTOCOL)

def parse_vcf(filename):
	#do the same as above for variant call file
	print "Starting to make vcf dictionary..."
	vcf = {}
	vcf2 = {}
	var_count = 0

	with open(filename) as f:
		for line in f:
			if line[:1] == "#":
				continue
			linedata = line.strip().split("\t")
			if linedata[6] == "PASS": #quality filter
				if linedata[0] not in vcf:
					vcf[linedata[0]] = []
				vcf[linedata[0]].append(int(linedata[1]))
				var_count += 1

	print var_count, "total variants in Mspr vcf."

	#same nested dictionaries as above
	for chrom in vcf:
		varlist = vcf[chrom]
		varlist.sort()
		splits = len(varlist)/100
		count = 0

		for i in range(0,len(varlist),splits):
			newlist = varlist[i:i+splits]
			maxx = max(newlist)
			minn = min(newlist)
			key = (minn,maxx)
			if chrom not in vcf2:
				vcf2[chrom] = {}
			vcf2[chrom][key] = {}
			
			newsplits = len(newlist)/50
			if newsplits == 0:
				vcf2[chrom][key][key] = newlist
			else:
				for j in range(0,len(newlist),newsplits):
					innerlist = newlist[j:j+newsplits]
					nmaxx = max(innerlist)
					nminn = min(innerlist)
					nkey = (nminn,nmaxx)
					vcf2[chrom][key][nkey] = innerlist


	for chrom in vcf2:
		newdict = {chrom:vcf2[chrom]}
		with open("/work/tkang/primary_ASE/vcf/mba/no_repeatmasker/bychrom/$"+chrom+"$vcf_bychrom.pkl","wb") as f:
			pickle.dump(newdict,f,pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
	
	#skip these chromosomes
	passed = ["X","Y","MT"]

	which = argv[1]

	if which == "gtrd":
		make_TF_dict("/bigrock_home/tkang/primary_ASE/gtrd/Mus_musculus_meta_clusters.interval")

	elif which == "vcf":
		parse_vcf("/work/tkang/primary_ASE/vcf/mba/no_repeatmasker/STFtoPWK_global_filtsnps.vcf")

	else:
		print "gtrd or vcf?"
