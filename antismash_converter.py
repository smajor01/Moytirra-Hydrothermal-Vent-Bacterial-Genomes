#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###########################################################################################################################################################################
#### This Python script will parse antiSMASH outputs in the directory where is executed.
#### It will transform index.html files into csv files
#### It will also produce a summary table with the counts of the BGCs types per genome: "bgc_type_table.csv"
#### Usage: python antismash_html_tocsv.py /path_with_antismash_files ########################################################################
###########################################################################################################################################################################

import argparse
import sys
parser=argparse.ArgumentParser(
	description='''This Python script searchs for antiSMASH html outputs inside 
	the indicated directory and produces to the Output directory 
	./antismash_Conversion_files :
	-csv file with antismash output for each genome: *_table.csv
	-Table with the information from all genomes concatenated: all_BGC_info.csv
	-Table with counts for the antismash scheme: bgc_antismash_class_counts.csv
	-Table with counts for the new classification scheme: BGCs_resumed.csv
	
	Note: necessary packages: simplified_scrapy, pandas
	pip install simplified-scrapy	
	''')

__author__ = 'Sandra Godinho Silva (sandragodinhosilva@gmail.com)'
__version__ = '0.5'
__date__ = '29-03-2022'

parser.add_argument('inputDirectory', help='Path to the input directory.')

###############################################################################
# Classification scheme: 
# NRPS: only nrps
NRPS = ["NRPS-like", "NRPS"]
# NRPS_other: nrps and other type, except pks
NRPS_other = ["thioamide-NRP", "NRPS,siderophore", "NRPS,lanthipeptide", 
				"terpene,NRPS-like,betalactone", "NRPS,indole", 
				"ladderane,NRPS","NRPS-like,bacteriocin","NRPS,proteusin,LAP"
				"NRPS-like,lanthipeptide", "NRPS,betalactone", 
				"NRPS-like,siderophore","NRPS,ladderane","NRPS-like,terpene", 
				"NRPS,LAP,proteusin", "arylpolyene,resorcinol,NRPS",
				"NRPS-like,betalactone", "NRPS,terpene", "siderophore,NRPS", 
				"terpene,NRPS-like", "NRPS-like,NRPS,siderophore",
				"NRPS-like,lanthipeptide", "NRPS,proteusin,LAP"
				]
				
# NRPS_PKS_hybrid: nrps with pks
NRPS_PKS_hybrid =["T1PKS,NRPS", "NRPS,T1PKS", "NRPS,transAT-PKS", 
				  "PKS-like,transAT-PKS,NRPS", "NRPS-like,T3PKS", 
				  "NRPS-like,T1PKS","NRPS,T1PKS,lanthipeptide",
				  "T3PKS,NRPS","NRPS,T1PKS,T3PKS", "T1PKS,NRPS-like", 
				  "transAT-PKS,transAT-PKS-like,NRPS-like,PKS-like,T3PKS", 
				  "T3PKS,NRPS,T1PKS", "T1PKS,NRPS", "NRPS,T1PKS,betalactone", 
				  "NRPS,bacteriocin","T3PKS,hglE-KS,siderophore,NRPS,T1PKS", 
				  "hglE-KS,T1PKS,NRPS,betalactone","transAT-PKS,NRPS",
				  "transAT-PKS-like,transAT-PKS,PKS-like,NRPS,T1PKS", 
				  "betalactone,NRPS-like", "transAT-PKS,NRPS-like", 
				  "NRPS,T1PKS,bacteriocin","NRPS-like,hglE-KS,T1PKS",
				  "NRPS,T1PKS,siderophore,hglE-KS,T3PKS", "T3PKS,NRPS-like",
				  "NRPS,T3PKS","NRPS,T1PKS,T3PKS","thioamide-NRP",
				  "transAT-PKS-like,transAT-PKS,T3PKS,PKS-like,NRPS-like",
				  "NRPS,hglE-KS,T1PKS", "transAT-PKS,NRPS,PKS-like",
				  "transAT-PKS,transAT-PKS-like,NRPS-like,PKS-like,T3PKS",
				  "T1PKS,hglE-KS,NRPS,siderophore", "transAT-PKS,NRPS,PKS-like",
				  "NRPS,T1PKS,siderophore","NRPS-like,T1PKS,NRPS"
				  ]
# transAT_PKS: only transAT_PKS (may have other types of PKS)
transAT_PKS = ["transAT-PKS","transatpks","transAT-PKS,PKS-like",
			   "transAT-PKS,PKS-like", 
			  "transAT-PKS5", "transAT-PKS-like", 
			  "transAT-PKS-like,transAT-PKS,PKS-like", 
			  "transAT-PKS,PKS-like,transAT-PKS-like",
			  "transAT-PKS,transAT-PKS-like", "transAT-PKS,bacteriocin",
			  "transAT-PKS-like,transAT-PKS,PKS-like,ladderane",
			  "transAT-PKS-like,transAT-PKS",
			  "transAT-PKS,PKS-like,ladderane"     
			  ]
# PKSI
PKSI = ["t1pks", "T1PKS"]
# PKSII
PKSII = ["t2pks"]
# PKSIII
PKSIII = [ "t3pks", "3PKS", "T3PKS", "T3PKS,betalactone", "T3PKS,arylpolyene", 
		   "arylpolyene,T3PKS", "arylpolyene,resorcinol,T3PKS", 
		   "lanthipeptide,T3PKS,bacteriocin",
		   "betalactone,T3PKS", "terpene,T3PKS", "T3PKS,terpene", 
		   "T3PKS,arylpolyene,resorcinol", "T3PKS,resorcinol"
		   ]
#PKS_other: combination of pks with other pks or with other types (except nrps)
PKS_other = ["otherks", "hglks", "PKS", "PKS-like", "hglE-KS", "hglE-KS,T1PKS",
			   "T1PKS,hglE-KS", "hglE-KS,T1PKS,terpene", 
			   "ladderane,transAT-PKS,PKS-like,transAT-PKS-like", "T1PKS,PUFA",
			   "ladderane,transAT-PKS,PKS-like","T1PKS,PUFA,hglE-KS"
			   "transAT-PKS,PKS-like,ladderane", "T1PKS,hglE-KS,terpene",
			   "lanthipeptide,T1PKS,hglE-KS", "hglE-KS,T1PKS,lanthipeptide",
			   "terpene,hglE-KS,PUFA,T1PKS", "terpene,T1PKS,hglE-KS",
			   "terpene,hglE-KS,T1PKS", "hglE-KS,PUFA,T1PKS", 
			   "PUFA,T1PKS,hglE-KS","hglE-KS,terpene,T1PKS","PUFA,T1PKS",
			   "lanthipeptide,hglE-KS,T1PKS", "hglE-KS,T1PKS,PUFA", 
			   ]
# Saccharides
Saccharides=["amglyccycl", "oligosaccharide", "cf_saccharide", "saccharide"]
# siderophore
Siderophore = ["siderophore"]
#Terpene
Terpene=["terpene"]
# only RiPPs
RiPPs= ["RiPP-like","RRE-containing","RRE-containing,RiPP-like","lantipeptide", "thiopeptide", "bacteriocin", "linaridin", "proteusin", 
		   "cyanobactin", "glycocin", "LAP", "lassopeptide", "sactipeptide", 
		   "bottromycin", "head_to_tail", "microcin", "microviridin", 
		   "lanthipeptide", "lipolanthine", "RaS-RiPP", "fungal-RiPP",
		   "bacteriocin,lanthipeptide", "lanthipeptide,bacteriocin",
		   "thiopeptide,LAP", "LAP,proteusin", "proteusin,LAP","RaS-RiPP",
		   "proteusin,LAP,bacteriocin","LAP,proteusin,bacteriocin",
		   "TfuA-related"
		   ]
Arylpolyene = ["arylpolyene"]

Resorcinol = ["resorcinol"]

# Others: diversified combinations and bgcs that don't fit previous classes
Others = ["acyl_amino_acids","aminocoumarin", "ectoine", 
			"butyrolactone", "nucleoside", "melanin", "phosphoglycolipid", 
			"phenazine", "phosphonate", "other", "cf_putative", 
			"indole", "ladderane", "PUFA", "furan", "hserlactone", "fused", 
			"cf_fatty_acid",  "blactam", "fatty_acid" "PpyS-KS", 
			"CDPS", "betalactone", "PBD", "tropodithietic-acid", "NAGGN", 
			"halogenated",  "terpene,bacteriocin","arylpolyene,bacteriocin",
			"arylpolyene,resorcinol", "resorcinol,arylpolyene", 
			"siderophore,terpene","terpene,ladderane","bacteriocin,acyl_amino_acids",
			"lanthipeptide,terpene", "arylpolyene,lanthipeptide,resorcinol",
			"acyl_amino_acids,bacteriocin", "ladderane,terpene", 
			"arylpolyene,resorcinol","arylpolyene,resorcinol,bacteriocin", 
			"lanthipeptide,siderophore","bacteriocin,arylpolyene,resorcinol", 
			"siderophore,bacteriocin", "terpene,lanthipeptide",
			"terpene,siderophore","terpene,arylpolyene,resorcinol",
			"terpene,betalactone", "bacteriocin,siderophore", 
			"terpene,bacteriocin,siderophore","terpene,arylpolyene"
			]
###############################################################################
# Import necessary Python modules
import os
from simplified_scrapy import SimplifiedDoc,utils
import pandas as pd

###############################################################################
# Import input folder
inputDirectory = sys.argv[1]
os.chdir(inputDirectory)
rootdir = os.getcwd()

# Output folder creation
output_dir = os.path.join(rootdir,"antismash_Conversion_files")

try:
	os.mkdir(output_dir)
except:
	pass
print("Output folder: " + str(output_dir))
print("")
###############################################################################
# Step 1: transform html files into csv ("_table.csv) for each genome

gbk_dir = os.path.join(output_dir,"Gbk_tables")

try:
	os.mkdir(gbk_dir)
except:
	pass

d = {}
for subdir, dirs, files in os.walk(rootdir):
	if "index.html" in files: #to confirm that is an antismash folder
		d_files = []
		for file in files:
			if ("gbk" in file) & ("out" not in file): #to get tbk files
				d_files.append(file.replace(".gbk",""))
			elif "index.html" in file: #to get html file
				print(file)
				name = str(subdir).split("/")[-1]
				index_path = os.path.join(subdir, file)
				html = utils.getFileContent(index_path) # Get data from file
				doc = SimplifiedDoc(html)
				rows = []
				tables = doc.selects('table.region-table')
				for table in tables:
					trs = table.tbody.trs
					for tr in trs:
						rows.append([td.text for td in tr.tds])
		d_files = sorted(d_files, key=str.lower)
		df = pd.DataFrame(rows)
		if df.empty:
			print("Genome without annotation")
			df.to_csv(os.path.join(gbk_dir, name + '_table.csv'), index=False)
		else:
			try:
				df.columns= ["Region","Type","From", "To","Known", "Known_type", "Known_similarity"]
			except:
				df.columns= ["Region","Type","From", "To","Known"]
			df_names = pd.DataFrame(d_files, columns=['Correct_name'])
			df = pd.concat([df, df_names], axis=1, sort=False)
			df.drop_duplicates(subset=["Region"], keep="first", inplace=True)
			df.to_csv(os.path.join(gbk_dir, name + '_table.csv'), index=False)
			print("Created: "+ name + '_table.csv')

###############################################################################                
# Step 2: Extract information from each genome _table file
print("STEP 2")
record_genomes_used = []
d_types = {}
d_bgc = {}               
df_major= pd.DataFrame()
for subdir, dirs, files in os.walk(gbk_dir):
	for file in files:
		if "_table.csv" in file:
			file_path = os.path.join(subdir, file)
			name = str(file).replace("_table.csv", "")
			record_genomes_used.append(name)
			print(file_path)
			try:
				df = pd.read_csv(file_path)
				if str(df.iloc[0,0]).startswith("Region"): # make sure the correct file is being parsed
					d_types[name] = df["Type"].tolist()
					df = df.assign(Genome=name)[['Genome'] + df.columns.tolist()]
					df_major = pd.concat([df_major, df], sort=False)
			except:
				print('csv file is empty for genome: ' + str(name))
###############################################################################
# Step3: Organize csv file that has all information
df_major.columns = ["Genome","Region_Type","antiSMASH_classif","From","To",
		"Most_similar_known_ cluster","Most_similar_classif","Similarity", 
		"Correct_name"]           

df_major["To"] = df_major["To"].str.replace(",","")
df_major["To"] = pd.to_numeric(df_major["To"])
df_major["From"] = df_major["From"].str.replace(",","")
df_major["From"] = df_major["From"].fillna("1")
df_major["From"] = pd.to_numeric(df_major["From"])

df_major["Size(bp)"] = df_major["To"] - df_major["From"]

cols = ["Genome","Region_Type","antiSMASH_classif","From","To","Size(bp)",
		"Most_similar_known_ cluster","Most_similar_classif",
		"Similarity", "new_classif", "Correct_name"] 

df_major = df_major.reindex(columns = cols)     

x=0
for x in (range(len(df_major))):
	a = df_major.iloc[x,2]
	if a in NRPS:
		df_major.iloc[x,9] = "NRPS"
		#print(a +" is NRPS")
	elif a in NRPS_other:
		df_major.iloc[x,9] ="NRPS_other"
		#print(a +" is NRPS_other")
	elif a in NRPS_PKS_hybrid:
		df_major.iloc[x,9] = "NRPS_PKS_hybrid"
		#print(a +" is NRPS_hybrid")
	elif a in transAT_PKS:
		df_major.iloc[x,9] = "transAT_PSK"
		#print(a +" is transAT_PKS")    
	elif a in PKSI:
		df_major.iloc[x,9] ="PKSI"
		#print(a +" is PKSI")
	elif a in PKSII:
		df_major.iloc[x,9] ="PKSII"
		#print(a +" is PKSI")
	elif a in PKSIII:
		df_major.iloc[x,9] ="PKSIII"
		#print(a +" is PKSIII")
	elif a in PKS_other:
		df_major.iloc[x,9] ="PKS_other"
		#print(a +" is PKS_other")
	elif a in Saccharides:
		df_major.iloc[x,9] ="Saccharides"
		#print(a +" is Saccharides")        
	elif a in Terpene:
		df_major.iloc[x,9]="Terpene"
		#print(a +" is terpene")
	elif a in Siderophore:
		df_major.iloc[x,9] ="Siderophore"
		#print(a +" is siderophore")            
	elif a in RiPPs:
		df_major.iloc[x,9] ="RiPPs"
	elif a in Arylpolyene:
		df_major.iloc[x,9] ="Arylpolyene"
	elif a in Resorcinol:
		df_major.iloc[x,9] ="Resorcinol"
	elif a in Others:
		df_major.iloc[x,9] ="Others"
	else:
		df_major.iloc[x,9] = "Unclassified"
		print(str(a) + " is unclassified")
	x +=1

cols = ["Genome","Region_Type","antiSMASH_classif", "new_classif","From","To",
		"Size(bp)","Most_similar_known_ cluster","Most_similar_classif",
		"Similarity", "Correct_name"]
df_major = df_major.reindex(columns = cols)    

df_major.to_csv(os.path.join(output_dir,"all_BGC_info.csv"), index=False) 
print("Table with all information concatenated was created: all_BGC_info.csv")
print("")
###############################################################################
# Step4: Create count table for antismash classification
from collections import Counter
df_types = pd.DataFrame({k:Counter(v) for k, v in d_types.items()}).fillna(0).astype(int)

df_types.to_csv(os.path.join(output_dir,r"bgc_antismash_class_counts.csv"))
print("Table with counts with the antismash scheme was created: bgc_antismash_class_counts.csv.csv")
print("It includes the antismash results from genomes: " + str(record_genomes_used))
print("")
###############################################################################
# Step5: Create count table for new classification
df_bgc_resumed = df_types.copy()
df_bgc_resumed = df_bgc_resumed.reset_index()

def BGC_classifier(df):
	x=0
	dd = {}
	for x in (range(len(df_bgc_resumed))):
		a = df_bgc_resumed.iloc[x,0]
		x+=1
		if a in NRPS:
			dd[a] = "NRPS"
		#print(a +" is NRPS")
		elif a in NRPS_other:
			dd[a] ="NRPS_other"
		#print(a +" is NRPS_other")
		elif a in NRPS_PKS_hybrid:
			dd[a] ="NRPS_PKS_hybrid"
		#print(a +" is NRPS_hybrid")
		elif a in transAT_PKS:
			dd[a] ="transAT_PSK"
		#print(a +" is transAT_PKS")    
		elif a in PKSI:
			dd[a] ="PKSI"
		#print(a +" is PKSI")
		elif a in PKSII:
			dd[a] ="PKSII"
		#print(a +" is PKSI")
		elif a in PKSIII:
			dd[a] ="PKSIII"
		#print(a +" is PKSIII")
		elif a in PKS_other:
			dd[a] ="PKS_other"
		#print(a +" is PKS_other")
		elif a in Saccharides:
			dd[a] ="Saccharides"
		#print(a +" is Saccharides")        
		elif a in Terpene:
		   dd[a] ="Terpene"
		#print(a +" is terpene")
		elif a in Siderophore:
			dd[a] ="Siderophore"
   		#print(a +" is siderophore")
		elif a in Arylpolyene:
			dd[a]  ="Arylpolyene"
		elif a in Resorcinol:
			dd[a]  ="Resorcinol"            
		elif a in RiPPs:
			dd[a] ="RiPPs"
		#print(a +" is RiPPNRPS = ["NRPS-like", "NRPS"]
		elif a in Others:
			dd[a] ="Others"
		else:
			dd[a] = "Unclassified"
			print(str(a) + " is unclassified")
	return dd

dd = BGC_classifier(df_bgc_resumed)        
df_bgc_resumed_t = df_bgc_resumed.set_index("index")
df_bgc_resumed_t.rename(index={k: v for k, v in dd.items()}, inplace=True)   
df_bgc_resumed_t = df_bgc_resumed_t.reset_index()

final = df_bgc_resumed_t.groupby("index", as_index=False).sum()
final.to_csv(os.path.join(output_dir,"BGCs_resumed.csv"), index=False)
print("Table with counts for the new classification scheme was created: BGCs_resumed.csv")

## END
