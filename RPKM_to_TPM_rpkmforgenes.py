#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  RPKM_to_TPM_stuartOutFrmt.py
#  
#  Copyright 2018 oscar <oscar@oscar-J53kOiSr>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import os

from numpy import all,array,divide,float32,ma,mean,multiply,sum,where,zeros

#########################
#########################
# Specific file
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("file", help="the file to be processed",
type=str)
args = parser.parse_args()

#print the square of user input from cmd line.
outFlRpkmforgenes = args.file
print 'Computing TPM for file %s...'%outFlRpkmforgenes


########
#Methods
########
def RPKM_GtoTPM_G_rpkmforgenes(ar_RPKM_G,R,r_l):
	"""
	Function to convert RPKM to TPM for a gene, specifically to account
	for density in rpkmforgenes.
	Input: ar_RPKM_G (reads per kilobase per million mapped reads) is an
	array with the RPKM values calculated for all genes in a library. R 
	is the total number of reads mapped (in the library). r_l is the 
	average read length.
	Output: ar_TPM_G (transcripts per million) is an array of TPMs 
	calculated for genes in the same order as the input ar_RPKM_G array. 
	"""
	sclngFctr = divide(R,1e9,dtype=float32)
	ar_D_G = multiply(ar_RPKM_G,sclngFctr)
	T_G = sum(multiply(ar_D_G,r_l))
	sclngFctr2 = multiply(r_l,1e6)
	nomntr = multiply(ar_D_G,sclngFctr2)
	TPM_G = divide(nomntr,T_G,dtype=float32)
	return TPM_G
	
def wrprRun_rpkmforgenes(outFlRpkmforgenes,outFlTpmforgenes,outFrmt='replace', \
	r_l=43):
	"""
	Function to run rpkm for genes and return TPMs
	Input: outFlRpkmforgenes is the output file with the results for 
	rpkmforgenes. outFlTpmforgenes is the output file with the results 
	for TPMs calculated for rpkmforgenes output. pthRpkforgenes is the 
	path to pthRpkforgenes. outFrmt is the format to output the 
	calculated values; these can be: "add": to just append the TPM value
	to the output file, "table_tpm" to just report the tpms (in 
	addition to the annotation), "replace" to replace rpkm for tpm.
	r_l is the average read length.
	"""
	#Assert output format is valid
	try:
		assert outFrmt in {'add','table_tpm','replace'}
	except:
		raise Exception('%s is not a valid output format'%outFrmt)
	#Test paths for input files
	crrntPth = os.getcwd()
	if not os.path.split(outFlRpkmforgenes)[0]:
		outFlRpkmforgenes = os.path.join(crrntPth,outFlRpkmforgenes)
		if os.path.exists(outFlRpkmforgenes):
			print Warning('%s is going to be overwritten'% \
			outFlRpkmforgenes)
	#assert input files exists
	try:
		assert os.path.exists(outFlRpkmforgenes)
	except:
		raise Exception('%s does not exist'%outFlRpkmforgenes)
	#Return data from output file
	lGnNms,lTrnscptCds,aSmpl_RPKM_G,aSmplCnts,aSmplName, \
	aSmplAllReads,aSmplMpdReads,argmnts = \
	prsOutFl_rpkmforgenes(outFlRpkmforgenes)
	#Filter TPMs for non-zeros
	aSmpl_RPKM_GNonZeros = ma.masked_less_equal(aSmpl_RPKM_G,0)
	aSmpl_TPM_G = zeros(aSmpl_RPKM_GNonZeros.shape,dtype=float32)
	#Compute TPMs
	for smplPos in xrange(len(aSmplName)):
		ar_RPKM_G = aSmpl_RPKM_G[smplPos]
		ar_RPKM_GNonZeros = aSmpl_RPKM_GNonZeros[smplPos]
		mpdReads = aSmplMpdReads[smplPos]
		TPM_GNonZeros = \
		RPKM_GtoTPM_G_rpkmforgenes(ar_RPKM_GNonZeros.compressed(), \
		mpdReads,r_l)
		posToFill = where(ar_RPKM_G>0)
		TPM_G = zeros(ar_RPKM_G.shape,dtype=float32)
		TPM_G[posToFill] = TPM_GNonZeros
		#Print ratio between ar_RPKM_G and TPM_G
		print 'sample,mean(RPKM/TPMs)->', \
		aSmplName[smplPos], \
		mean(ar_RPKM_GNonZeros.compressed()/TPM_GNonZeros)
		aSmpl_TPM_G[smplPos] = TPM_G
	#Write output
	outFlopnd = open(outFlTpmforgenes,'w')
	if outFrmt in {'add','replace'}:
		outFlopnd.write \
		('#samples\t%s\n#allmappedreads\t%s\n#normalizationreads\t%s\n#arguments\t%s\n' \
		%('\t'.join(aSmplName),'\t'.join(aSmplAllReads), \
		'\t'.join(str(v) for v in aSmplMpdReads),argmnts))
	for gnPos in xrange(len(lGnNms)):
		lOut = [lGnNms[gnPos],lTrnscptCds[gnPos]]
		if outFrmt in {'add'}:
			for smplPos in xrange(len(aSmplName)):
				lOut.append(str(aSmpl_RPKM_G[smplPos][gnPos]))
			for smplPos in xrange(len(aSmplName)):
				lOut.append(str(aSmplCnts[smplPos][gnPos]))
			for smplPos in xrange(len(aSmplName)):
				lOut.append(str(aSmpl_TPM_G[smplPos][gnPos]))
		elif outFrmt in {'replace'}:
			for smplPos in xrange(len(aSmplName)):
				lOut.append(str(aSmpl_TPM_G[smplPos][gnPos]))
			for smplPos in xrange(len(aSmplName)):
				lOut.append(str(aSmplCnts[smplPos][gnPos]))
		else:#table_tpm
			for smplPos in xrange(len(aSmplName)):
				lOut.append(str(aSmpl_TPM_G[smplPos][gnPos]))
		outFlopnd.write('%s\n'%'\t'.join(lOut))
	outFlopnd.close()
	return 0
	
def prsOutFl_rpkmforgenes(outFlRpkmforgenes):
	"""
	Function to parse the output files generated with rpkmforgenes
	Input: outFlRpkmforgenes is an output file from rpkmforgenes that 
	includes the name of the sample, the reads mapped to coding regions, 
	and a list of gene names, ncbi transcript codes, rpkm and counts.
	Output: The method obtains the number of reads mapped to coding 
	regions from the input file, 
	"""
	aSmplName,aSmplAllReads,aSmplMpdReads,argmnts = None,None,None,None
	lGnNms,lTrnscptCds,aSmpl_RPKM_G,aSmplCnts = [],[],[],[]
	for l in open(outFlRpkmforgenes,'r'):
		if l.strip():
			if l[0]=='#':
				if l.split()[0]=='#samples':
					aSmplName = array([s for s in \
					l.splitlines()[0].split('\t')[1:] if s.strip()])
				elif l.split()[0]=='#allmappedreads':
					aSmplAllReads = array([a for a in \
					l.splitlines()[0].split('\t')[1:] if s.strip()])
				elif l.split()[0]=='#normalizationreads':
					aSmplMpdReads = array([float32(m) for m in \
					l.splitlines()[0].split('\t')[1:] if s.strip()])
				elif l.split()[0]=='#arguments':
					argmnts = l.splitlines()[0].split()[1]
			else:
				l=l.splitlines()[0].split('\t')
				lGnNms.append(l[0])
				lTrnscptCds.append(l[1])
				aSmpl_RPKM_G.append([float32(v) for v in \
				l[2:len(aSmplName)+2]])
				aSmplCnts.append([v for v in \
				l[len(aSmplName)+2:]])
	try:
		assert aSmplName is not None
	except:
		raise Exception('No sample name found in rpkmforgenes output file')
	try:
		assert aSmplAllReads is not None
	except:
		raise Exception('No total reads found in rpkmforgenes output file')
	try:
		assert aSmplMpdReads is not None
	except:
		raise Exception('No mapped reads found in rpkmforgenes output file')
	if len(lGnNms)==0:
		raise Exception('Gene expression is empty in the output file')
	lGnNms,lTrnscptCds,aSmpl_RPKM_G,aSmplCnts = array(lGnNms), \
	array(lTrnscptCds),array(aSmpl_RPKM_G),array(aSmplCnts)
	aSmpl_RPKM_G = aSmpl_RPKM_G.T
	aSmplCnts = aSmplCnts.T
	return lGnNms,lTrnscptCds,aSmpl_RPKM_G,aSmplCnts,aSmplName, \
	aSmplAllReads,aSmplMpdReads,argmnts



########
#Input and run 
########
#set input and output
#~ outFlRpkmforgenes = 'ensembl_rpkms_DEC13A.txt'#name of input file from (output from rpkmforgenes)
#~ outFlTpmforgenes = 'ensembl_rpkms_DEC13A.tpm'#name of output file with tpm values instead of rpkms (calculated as Wagner and Lynch 2012)
outFlTpmforgenes = '%s.tpm'%outFlRpkmforgenes#name of output file with tpm values instead of rpkms (calculated as Wagner and Lynch 2012)
#run
wrprRun_rpkmforgenes(outFlRpkmforgenes,outFlTpmforgenes)

#Log out for AUG13A
#sample,mean(RPKM/TPMs)-> SCG_WT_1 0.5111476
#sample,mean(RPKM/TPMs)-> SCG_WT_2 0.5406556
#sample,mean(RPKM/TPMs)-> SCG_WT_3 0.52702093
#sample,mean(RPKM/TPMs)-> SCG_EGLN3_KO_1 0.54187053
#sample,mean(RPKM/TPMs)-> SCG_EGLN3_KO_2 0.5245067
#sample,mean(RPKM/TPMs)-> SCG_EGLN3_KO_3 0.5623313

#Log out for DEC13A
#sample,mean(RPKM/TPMs)-> SCG_E3_KO_4 0.5313835
#sample,mean(RPKM/TPMs)-> SCG_E3_KO_5 0.5186138
#sample,mean(RPKM/TPMs)-> SCG_E3_WT_4 0.58979213
#sample,mean(RPKM/TPMs)-> SCG_E3_WT_5 0.5523013
#sample,mean(RPKM/TPMs)-> SCG_E3_WT_6 0.6113889
#sample,mean(RPKM/TPMs)-> SCG_E3_WT_7 0.5682907
