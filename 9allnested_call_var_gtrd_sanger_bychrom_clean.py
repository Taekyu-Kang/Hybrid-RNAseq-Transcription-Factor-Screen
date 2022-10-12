from sys import argv
import numpy as np
import time
import pickle
import os
import collections

#combines both outputs from gtrd_sanger_bychrom.py
#finds sequence variants located within reported transcription factor binding sites

def load_pkl(filename):
	with open(filename,"rb") as f:
		return pickle.load(f)

#order dictionary keys to speed up search
def order_dict(dat):
	newfinal = {dat.keys()[0]:""}
	dat2 = dat[dat.keys()[0]]
	newdat = collections.OrderedDict(sorted(dat2.items()))

	for i in newdat:
		newdat[i] = collections.OrderedDict(sorted(newdat[i].items()))

	newfinal[dat.keys()[0]] = newdat

	return newfinal

#find all transcription factor binding sites located in a 5kb window upstream of each gene's transcriptional start site
def count_TF(dat_list,TF):

	TFs = []
	chrom0 = dat_list[0][2:]
	sites = TF[chrom0]
	strand = dat_list[1]
	start = dat_list[2]
	end = dat_list[3]

	if strand == "+":

		for site in sites: #look for first window that gene TSS is in
			if site[1] < start:
				continue
			if start in range(site[0],site[1]):
				start_window = sites[site]
				break

		for site2 in start_window: #second window
			if site2[1] < start:
				continue
			if start in range(site2[0],site2[1]):
				sw2 = start_window[site2]
				break

		for pos in sw2: #find all TF binding sites in the upstream window
			if len(pos) == 0:
				continue
			if pos[0] in range(start-5000,start) or pos[1] in range(start-5000,start):
				TFs.append(pos)

	elif strand == "-": #same as above for reverse strand

		for site in sites:
			if site[0] > end:
				continue
			if end in range(site[0],site[1]):
				end_window = sites[site]
				break

		for site2 in end_window:
			if site2[0] > end:
				continue

			if end in range(site2[0],site2[1]):
				ew2 = end_window[site2]
				break

		for pos in ew2:
			if len(pos) == 0:
				continue
			if pos[0] in range(end,end+5000) or pos[1] in range(end,end+5000):
				TFs.append(pos)

	return TFs

#look for variants in the reported binding sites upstream of gene TSS
def check_vars(new_mm,gtrd,vcf):
	final_dict = {}
	count = 0
	varsinTF = 0
	
	print "Starting variant counting..."

	all_start = time.time()

	for gene in new_mm:

		if new_mm[gene][0] == "mmY" or new_mm[gene][0] == "mmX" or new_mm[gene][0] == "mmMT": #not interested in genes that fall in these regions
			continue

		#count += 1
		#if count % 5000 == 0:
		#	print count, "genes processed of", len(new_mm), "( took", (time.time()-all_start)/60, "min )"
		#can implement a counter to see how long the script will take, but this will slow it down

		if gene not in final_dict:
			final_dict[gene] = {}
		gene_TFs = count_TF(new_mm[gene],gtrd)
		
		instart = time.time()
		chrom = new_mm[gene][0]

		for site in gene_TFs: #go through each upstream TF site, look for variants that fall within these sites
			start = site[0]
			end = site[1]
			name = site[2]

			if name not in final_dict[gene]:
				final_dict[gene][name] = 0

				for var in vcf[chrom]: #find first window for variants
					if var[1] < start:
						continue
					if start in range(var[0],var[1]):
						start_window = vcf[chrom][var]
						#print start_window
					if end in range(var[0],var[1]):
						end_window = vcf[chrom][var]
						break

				if start_window == end_window: #if the binding site does not span two windows, only need to search one set of nested windows
					for var2 in start_window:
						if var2[1] < start:
							continue
						if start in range(var2[0],var2[1]):
							sw2 = start_window[var2]
						if end in range(var2[0],var2[1]):
							ew2 = start_window[var2]
							break

					if sw2 == ew2:
						for pos in sw2:	
							if pos in range(start,end+1):
								final_dict[gene][name] += 1
								varsinTF += 1

					else:
						aw2 = [sw2,ew2]
						for w in aw2:
							for pos in w:
								if pos in range(start,end+1):
									final_dict[gene][name] += 1
									varsinTF += 1

				else:
					for var2 in start_window:
						if var2[1] < start:
							continue
						if start in range(var2[0],var2[1]):
							sw2 = start_window[var2]
							break

					for var2 in end_window:
						if var2[1] < end:
							continue
						if end in range(var2[0],var2[1]):
							ew2 = end_window[var2]
							break

					aw2 = [sw2,ew2]
					for w in aw2:
						for pos in w:
							if pos in range(start,end+1):
								final_dict[gene][name] += 1
								varsinTF += 1
			
			else:
				continue

	print "Done counting variants in TF binding sites."
	print varsinTF, "total variants in TF binding sites of 38,589,139."

	with open("/work/tkang/primary_ASE/ctPEAK_GTRD_allTF_Sanger_varct.pkl","wb") as f:
		pickle.dump(final_dict,f,pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
	gtrd_dir = "/work/tkang/primary_ASE/gtrd/bychrom/"
	vcf_dir = "/work/tkang/primary_ASE/vcf/mba/no_repeatmasker/bychrom/"
	start1 = time.time()
	all_gtrd = {}
	all_vcf = {}

	for file in os.listdir(gtrd_dir):
		gtrd = load_pkl(gtrd_dir+file)
		newgtrd = order_dict(gtrd)
		all_gtrd[newgtrd.keys()[0]] = newgtrd[newgtrd.keys()[0]]

	print "Done loading gtrd files,", len(all_gtrd), "chromosomes,", all_gtrd.keys()
	print "Took", time.time()-start1, "sec."

	start2 = time.time()

	for file in os.listdir(vcf_dir):
		vcf = load_pkl(vcf_dir+file)
		newvcf = order_dict(vcf)
		all_vcf[newvcf.keys()[0]] = newvcf[newvcf.keys()[0]]


	print "Done loading vcf files,", len(all_vcf), "chromosomes,", all_vcf.keys()
	print "Took", time.time()-start2, "sec."

	annot = load_pkl("/bigrock_home/tkang/primary_ASE/Mmus_geneID_strand_dict.pkl")
	print "Done loading annotations."

	check_vars(annot,all_gtrd,all_vcf)

