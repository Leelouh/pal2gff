#!/usr/bin/env python3

#usage : pal2gff.pl palindrome_file.out > file.gff

#lancement de palindrome :
##palindrome -sequence seq.fa -minpallen 5 -maxpallen 20 -gaplimit 6 -nummismatches 0 -outfile truc.pal -nooverlap

import sys,re,os,cgi,csv, getopt, argparse, time, math, shutil
from copy import deepcopy
from random import randrange, uniform
from numpy import *
import scipy as sp
from pandas import *

def create_output_dir (output_dir) : 
	directory = output_dir
	if not os.path.isdir(directory):
		os.makedirs(directory)

def pal_2_bed (input_file) :
	lines = input_file.split("\n")
	i = 0

	for line in lines :
		search_ID = re.search("Palindromes\sof\:", line) #search id/chr for each event of palindrom
		if search_ID : 
			i = i + 1
			line_splitted = re.split('\s+', line)
			line_splitted = list(filter(None, line_splitted))
			chr = line_splitted[2]
			while (i<len(lines)) :
				search_pos = re.search('^\d', lines[i]) #search for position (start and end)
				#search_aln = re.search('\|', lines[j]) #for calcul of mismatch nb between each hits
				search_next_ID = re.search("Palindromes\sof\:", lines[i]) #stop if we have "palindrome of"
				if (search_next_ID):
					break
				elif (search_pos):
					search_aln = re.search('\|', lines[i+1])
					if (search_aln):
						hit_splitted_line1 = re.split('\s+', lines[i])
						hit_splitted_line2 = re.split('\s+', lines[i+2])
						pos_start = int(hit_splitted_line1[0])
						pos1_spacer = int(hit_splitted_line1[2])
						pos2_spacer = int(hit_splitted_line2[2])
						pos_end = int(hit_splitted_line2[0])
						R1 = len(hit_splitted_line1[1])
						R2 = len(hit_splitted_line2[1])
						hit = re.sub(r"\s+", "", lines[i+1])

						hit_nb = len(hit)
						spacer = pos2_spacer - pos1_spacer - 1

						mis = R1 - hit_nb
						if (R1 == 6 and mis ==1):
							pass
						elif (R1 == 6 and mis == 0 and spacer in range(0,4)) :
							pass
						elif (R1 == 7 and mis == 0 and spacer in range(0,2)) :
							pass
						elif (R1 == 7 and mis == 1 and spacer in range(0,4)) :
							pass
						else :
							if (R1 == 8 and spacer in range(0,2)):
								goodAln = re.search('^\s*\|\|\|\|\|\|', lines[i+1])
								if (not goodAln):
									pass
								else :
									name_list=[str(chr),str(pos_start),str(pos_end)]
									name = '_'.join(name_list)
									ID = "ID="+name+";spacer="+str(spacer)+";repeat="+str(R1)+";mismatch="+str(mis)
									print (chr,"pal2gff","IR",pos_start,pos_end,".","+",".",str(ID),sep='\t')
							else :
								name_list=[str(chr),str(pos_start),str(pos_end)]
								name = '_'.join(name_list)
								ID = "ID="+name+";spacer="+str(spacer)+";repeat="+str(R1)+";mismatch="+str(mis)
								print (chr,"pal2gff","IR",pos_start,pos_end,".","+",".",str(ID),sep='\t')
				i = i+1

with open(sys.argv[1], 'r') as input_file_data : 
	input_file = input_file_data.read()

#create_output_dir(sys.argv[2]) 

pal_2_bed(input_file) 
