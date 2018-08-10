#!/usr/bin/env python
# coding=utf-8
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
from Bio import SeqIO
import networkx as nx

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

g=nx.MultiDiGraph()

## Adding all vertexes to the graph
with open(args.fasta, "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        g.add_node(record.id,
        seq=str(record.seq),
        length=len(record.seq))
handle.close()

# Adding edges from Ragout :
scaffolding_type_1="synteny"
current_scaffold=""
fromName=""
gap_length=0
with open(args.phylo) as f:
    for line in f:
        if line[0] != "#":
            line.rstrip("\n")
            elems=line.split("\t")
            if current_scaffold == elems[0]:
                if elems[4] == "N":
                    gap_length=elems[5]
                else:
                    g.add_edge(fromName,elems[5],type=scaffolding_type_1,fromStrand=from_orientation,ToStrand=elems[8],gap=gap_length, readsCount=-1)
                    fromName=elems[5]
                    from_orientation=elems[8]
            else:
                current_scaffold=elems[0]
                fromName=elems[5]
                from_orientation=elems[8]
f.close()

scaffolding_type_2="rna-seq"
with open(args.rna) as f:
    for line in f:
        line.rstrip("\n")
        elems=line.split("\t")
        fromName=elems[1]
        toName=elems[4]
        from_orientation=elems[3]
        to_orientation=elems[6]
        # Here we add the edge, but try to simplify it if it already exist because of Ragout
        numFromTo=g.number_of_edges(fromName,toName)
        edgeToAdd=True
        if numFromTo>0:
            for i in range(0,numFromTo):
                # Orientation consistency ?
                if ((g.edges[fromName,toName,i]['fromStrand']==fromStrand) and (g.edges[fromName,toName,i]['ToStrand']==ToStrand)) or ((g.edges[toName,fromName,i]['fromStrand']==strand_reverse(fromStrand)) and (g.edges[toName,fromName,i]['ToStrand']==strand_reverse(ToStrand))):
                    if g.edges[fromName,toName,i]["type"]!= scaffolding_type_2:
                        g.edges[fromName,toName,i]['readsCount']+=1
                        edgeToAdd=False
                    else:
                        g.edges[fromName,toName,i]["type"]=scaffolding_type_1 + "_AND_" + scaffolding_type_2
                        g.edges[fromName,toName,i]['readsCount']=1
                        edgeToAdd=False
                if not edgeToAdd:
                    break
                    
            if edgeToAdd:
                g.add_edge(fromName,toName,type=scaffolding_type_2,fromStrand=fromStrand,ToStrand=ToStrand,gap=-1,reads=1)

f.close()

nx.write_gml(g, args.outDir+"/graph.gml")