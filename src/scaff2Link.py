#!/usr/bin/env python
import argparse
version='SEDMATCHGITVERSION'
year=2018
authors=['Martin Binet','Julien Fouret']
contact='julien.fouret@fouret.me'
parser = argparse.ArgumentParser(description="""

#-----------------------|Scaff2link|-----------------------#
     Scaffolding with a linkage graph using 2 sources
------------------------------------------------------------
This software allows the importing of 2 types of edges in a 
Scaffolding graph depending on the source of the link. Cur-
rently only phylogenetic links from ragout software and rna
-seq based links from agouti are supported.

""",epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+", ".join(authors)+"\nfor more informations or enquiries please contact :\n"+contact,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-fasta', metavar='/path', required=True, help="""
fasta file with contigs or solid-scaffold to analyze""")
parser.add_argument('-phylo', metavar='/path', required=True, help="""
linkage map from ragout (format to specify)""")
parser.add_argument('-rna', metavar='/path', required=True, help="""
linkage map from agouti (format to specify)""")
parser.add_argument('-tree', metavar='tree.nh', required=True, help="""
Tree file with newick format""")
parser.add_argument('-outDir', metavar='/path', required=True, help="""Output directory""")
parser.add_argument('-polarity_first',metavar='priority',choices=["phylo","rna","None"],default="None",required=False,help="""
what type of link keep if there in case of node bipolarity""")
parser.add_argument('-dag_first',metavar='priority',choices=["phylo","rna","None"],default="None",required=False,help="""
what type of link keep if there in case of a complexe dag 
(not all node in the longest path)""")
#what type of link keep if there in case of cycle""")
#parser.add_argument('-fork_first',metavar='priority',choices=["phylo","rna","None"],default="None",required=False,help="""
#what type of link keep if there in case of a fork""")
parser.add_argument('-v',action='store_true', help="verbose")
args=parser.parse_args()

import sys
import datetime

def logger(ltype,msg):
	st = datetime.datetime.now()
	if (ltype=="err" or ltype=="warn"):
		sys.stderr.write(" :: ".join([st,ltype,msg]))
		if ltype=="err":
			sys.exit(1)
	elif ltype=="out":
		sys.stdout.write(" :: ".join([st,ltype,msg]))
	else:
		logger("warn","error in the code the logger function is called with a base type...")

if args.v:
	logger("out","Begin the analysis with Scaff2link version "+version)


# pour avoir les arguments c'est args.fasta / args.phylo / args.rna etc ... ça te renvoie la string (type string)