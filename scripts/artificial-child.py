#!/usr/bin/env python3

import sys
from collections import defaultdict, Counter, namedtuple
import vcf
import random
import math
import datetime

RecombinationMapEntry = namedtuple('RecombinationMapEntry', ['position', 'cum_distance'])

def interpolate(point, start_pos, end_pos, start_value, end_value):
	assert start_pos <= point <= end_pos
	if start_pos == point == end_pos:
		assert start_value == end_value
		return start_value
	return start_value + ((point - start_pos) * (end_value - start_value) / (end_pos - start_pos))

def recombination_cost_map(genetic_map, positions):

	assert len(genetic_map) > 0


	# Step 1: compute cumulative genetic distances from start of chromosome 
	#         to each position.
	cumulative_distances = []
	# i and j are such that genetic_map[i].position <= position <= genetic_map[j].position
	# i and j are None if no such values exist (because we are at the end of the list)
	i = None
	j = 0

	for position in positions:
		# update i to meet the invariant
		if (i == None) and (genetic_map[0].position <= position):
			i = 0
		while (i != None) and (i+1 < len(genetic_map)) and (genetic_map[i+1].position <= position):
			i += 1

		# update j to meet the invariant
		while (j != None) and (genetic_map[j].position < position):
			if j+1 < len(genetic_map): 
				j += 1
			else:
				j = None

		# interpolate
		if i == None:
			assert j != None
			d = interpolate(position, 0, genetic_map[j].position, 0, genetic_map[j].cum_distance)
		elif j == None:
			# Point outside the genetic map --> extrapolating using average recombination rate
			avg_rate = genetic_map[-1].cum_distance / genetic_map[-1].position
			d = genetic_map[-1].cum_distance + (position - genetic_map[-1].position) * avg_rate
		else:
			assert genetic_map[i].position <= position <= genetic_map[j].position
			d = interpolate(position, genetic_map[i].position, genetic_map[j].position, genetic_map[i].cum_distance, genetic_map[j].cum_distance)
		cumulative_distances.append(d)

	# Step 2: compute costs (= phred-scaled recombination probabilities between two positions)
	#result = [0]
	#for i in range(1, len(cumulative_distances)):
	#	d = cumulative_distances[i] - cumulative_distances[i-1]
	#	result.append(round(centimorgen_to_phred(d)))

	return cumulative_distances




def load_genetic_map(filename):
	genetic_map = []

	with open(filename,'r') as fid:

		# read and ignore first line
		fid.readline()

		# for each line only store the first and third value in two seperate list
		for line in fid:
			line_spl = line.strip().split()
			assert len(line_spl) == 3
			genetic_map.append(
				RecombinationMapEntry(position=int(line_spl[0]), cum_distance=float(line_spl[2]))
			)
	
	return genetic_map


gm_filename = sys.argv[1]
vcf_mother_filename = sys.argv[2]
vcf_father_filename = sys.argv[3]
child_sample_name = sys.argv[4]
mother_recomb_output_filename = sys.argv[5]
father_recomb_output_filename = sys.argv[6]
#snp_positionsm = sys.argv[2]

#snp_positionsc = sys.argv[4]

#snp_positionsf = sys.argv[5]

#truebreaksm= sys.argv[6]
#truebreaksf= sys.argv[7]
gm = load_genetic_map(gm_filename)
positions=[]

vcf_reader = vcf.Reader(filename=vcf_father_filename)

for record in vcf_reader:    
	pos=int(record.POS)
	positions.append(int(pos))

cumgd=recombination_cost_map(gm,positions)

prob = [0]
for i in range(1, len(cumgd)):
	d = cumgd[i] - cumgd[i-1]
	p = (1.0-math.exp(-(2.0*d)/100))/2.0
	prob.append(p)

#Mother
alle_m=[]
vcf_reader = vcf.Reader(filename=vcf_mother_filename)
i=0
flag=0
#m1=open("mother_vcf.vcf","w")
mp=open(mother_recomb_output_filename,"w")
#print("##fileformat=VCFv4.1",file=m1)
#print("##fileDate=23122015_13h48m54s",file=m1)
#print("##source=SHAPEIT2.v837",file=m1)
#print("##log_file=shapeit_23122015_13h48m54s_e2ce506c-6a31-4d8c-8efc-f25b46265ce8.log",file=m1)
#print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased Genotype\">",file=m1)
#print("#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG004",file=m1)
for record in vcf_reader:
	gt = record.samples[0]['GT']
	allele1=gt.split('|')
	x=random.random()
	if(x <= prob[i]):
		print(record.POS,file=mp)
		if(flag==0):
			flag=1
		else:
			flag=0
	if(flag==1):
		temp=allele1[0]
		allele1[0]=allele1[1]
		allele1[1]=temp
	i=i+1
	#al=str(allele1[0])+"|"+str(allele1[1])
	alle_m.append(allele1[0])
	#print(record.CHROM,"	",record.POS,"	.	",record.REF,"	",record.ALT,"	",record.QUAL,"	PASS	",record.FILTER,"	GT	",al,file=m1)
	
#father
alle_f=[]
vcf_reader = vcf.Reader(filename=vcf_father_filename)
i=0
flag=0
#f1=open("father_vcf.vcf","w")
fp=open(father_recomb_output_filename,"w")
#print("##fileformat=VCFv4.1",file=f1)
#print("##fileDate=23122015_13h48m54s",file=f1)
#print("##source=SHAPEIT2.v837",file=f1)
#print("##log_file=shapeit_23122015_13h48m54s_e2ce506c-6a31-4d8c-8efc-f25b46265ce8.log",file=f1)
#print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased Genotype\">",file=f1)
#print("#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG003",file=f1)
for record in vcf_reader:
	gt = record.samples[0]['GT']
	allele1=gt.split('|')
	x=random.random()
	if(x <= prob[i]):
		print(record.POS,file=fp)
		if(flag==0):
			flag=1
		else:
			flag=0
	if(flag==1):
		temp=allele1[0]
		allele1[0]=allele1[1]
		allele1[1]=temp
	i=i+1
	#al=str(allele1[0])+"|"+str(allele1[1])
	alle_f.append(allele1[0])
	#print(record.CHROM,"	",record.POS,"	.	",record.REF,"	",record.ALT,"	",record.QUAL,"	PASS	",record.FILTER,"	GT	",al,file=f1)

vcf_reader = vcf.Reader(filename=vcf_father_filename)
i=0
print("##fileformat=VCFv4.1")
print('##fileDate={}'.format(datetime.datetime.now().strftime('%Y%m%d')))
print("##source=artificial-child.py")
print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased Genotype\">")
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(child_sample_name))
for record in vcf_reader:
	al=str(alle_m[i])+"|"+str(alle_f[i])
	assert len(record.ALT) == 1
	print(record.CHROM,record.POS,".",record.REF,record.ALT[0],".","PASS",".","GT",al, sep='\t')
	i=i+1

