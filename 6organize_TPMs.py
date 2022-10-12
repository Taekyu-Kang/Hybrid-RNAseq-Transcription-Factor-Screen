from sys import argv
import pickle
import time
import pandas as pd
from functools import reduce
import os

#meta data file, will be dictionary with file name as key, and annotations as values
#htseq out files

def parse_filetodict(filename,sep="\t"):
	fdict = {}
	with open(filename) as f:
		for line in f:
			linedata = line.strip().split(sep)
			fdict[linedata[0]] = linedata[1]
	return fdict

#pandas filter for gene counts
def parse_htseq(filename,metadata,sep="\t"):

	fileID = filename.split("/")[-1]
	sampleID = metadata[fileID]

	htseq = pd.read_csv(filename,sep=sep,header=None)
	htseq = htseq.rename(columns={0:"Gene",1:sampleID})
	htseq1 = htseq[[x.startswith("EN") or x.startswith("MG") for x in htseq["Gene"]]]

	return htseq1

#function to convert counts to TPM
def make_TPM(colname,df):
	df[colname+"_len"] = df[colname].div(df.Length,axis=0)
	df["ssum"] = df[colname+"_len"].sum()
	df[colname+"_div"] = df[colname+"_len"].div(df.ssum,axis=0)
	df[colname+"_TPM"] = df[colname+"_div"].multiply(1000000)

#function to parse file containing homologous genes 
def parse_homs(filename,sep="\t"):
	homsdict = {"Gene":[],"GeneID":[]}
	with open(filename) as f:
		for line in f:
			linedata = line.strip().split(sep)
			if linedata[5] == "ortholog_one2one":
				homsdict["Gene"].append(linedata[0])
				homsdict["Gene"].append(linedata[3])
				homsdict["GeneID"].append(linedata[1])
				homsdict["GeneID"].append(linedata[1])
	homsdf = pd.DataFrame(homsdict)
	return homsdf

if __name__ == "__main__":

	outdir = argv[1]
	metadata_file = argv[2]
	transcript_length_file = argv[3]
	homs_file = argv[4]
	htseq_dir = argv[5] #name of this directory will be the name of the output files

	outID = htseq_dir.split("/")[-2]

	htseq_files = []
	for i in os.listdir(htseq_dir):
		if i.endswith("htseq.out"):
			htseq_files.append(htseq_dir+i)

	metadata_dict = parse_filetodict(metadata_file)
	print len(metadata_dict), "files to be processed."

	homs_df = parse_homs(homs_file)
	print len(homs_df), "genes with homology between M.mus and M.spr."

	#read in file containing transcript lengths collected from GTF
	transcript_length = pd.read_csv(transcript_length_file,sep="\t")
	print len(transcript_length), "genes with transcript length in GTF."

	htseq_all = []

	for htseq_file in htseq_files:
		htseq_dat = parse_htseq(htseq_file,metadata_dict)
		htseq_all.append(htseq_dat)

	htseq_all.append(transcript_length)

	#merge all htseq files along with transcript length
	combined = reduce(lambda left,right: pd.merge(left,right,on=["Gene"],how="outer"),htseq_all).fillna(0)
	print len(htseq_all), "files combined."
	
	combined = combined.set_index("Gene")
	combined1 = combined[(combined != 0).all(axis=1)] #removes those with 0 in all samples
	col_names = combined1.columns.tolist()
	col_names.remove("Length")

	#convert to TPM
	for col_name in col_names:
		make_TPM(col_name,combined1)

	#just grab TPM data
	col_names = combined1.columns.tolist()
	tpm_df = combined1[[name for name in col_names if "TPM" in name]]

	tpm_df.to_csv(outdir+outID+"_TPM.txt",sep="\t",index=True)

	tpm_df_homs = pd.merge(homs_df,tpm_df,on=["Gene"],how="outer").fillna(0)

	#split up allele specific data, then combine according to homologous genes and clean up
	mm_df = tpm_df_homs[tpm_df_homs["Gene"].str.contains("ENS")]
	ms_df = tpm_df_homs[tpm_df_homs["Gene"].str.contains("MGP")]

	mm_df = mm_df.drop("Gene",axis=1)
	mm_df = mm_df.set_index("GeneID")
	mm_newcols = [name+"_mm" for name in mm_df.columns.tolist() if "TPM" in name]
	mm_df.columns = mm_newcols

	ms_df = ms_df.drop("Gene",axis=1)
	ms_df = ms_df.set_index("GeneID")
	ms_newcols = [name+"_ms" for name in ms_df.columns.tolist() if "TPM" in name]
	ms_df.columns = ms_newcols

	org_df = pd.merge(mm_df,ms_df,on=["GeneID"],how="outer").fillna(0)
	org_colnames = org_df.columns.tolist()
	org_colnames.sort()
	sorted_df = org_df[org_colnames]

	sorted_df.to_csv(outdir+outID+"_organized_TPM.txt",sep="\t",index=True)
	print "Done generating organized TPM dataset."
	
