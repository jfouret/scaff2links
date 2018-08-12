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
linkage map from ragout (.*_scaffolds.agp)""")
parser.add_argument('-rna', metavar='/path', required=True, help="""
joint pairs from agouti (.*join_pairs.noise_free.txt)""")
parser.add_argument('-min_reads', metavar='N', required=False, default="10" , help="""
minimum number of reads to consider a gene scaffold""")
parser.add_argument('-rna_gap', metavar='N', required=False, default="1000" , help="""
gap length to use chen there is only rna-seq support""")
parser.add_argument('-lib_type', metavar='fr|rf|ff|rr', choices=["fr","rr","ff",'rf'],default="fr",required=False, help="""
Orientation of pairs""")

parser.add_argument('-outDir', metavar='/path', required=True, help="""Output directory""")



#parser.add_argument('-polarity_first',metavar='priority',choices=["phylo","rna","None"],default="None",required=False,help="""
#what type of link keep if there in case of node bipolarity""")
#parser.add_argument('-dag_first',metavar='priority',choices=["phylo","rna","None"],default="None",required=False,help="""
#what type of link keep if there in case of a complexe dag 
#(not all node in the longest path)""")

#what type of link keep if there in case of cycle""")
#parser.add_argument('-fork_first',metavar='priority',choices=["phylo","rna","None"],default="None",required=False,help="""
#what type of link keep if there in case of a fork""")
parser.add_argument('-v',action='store_true', help="verbose")
parser.add_argument('-d',action='store_true', help="Debug info")
args=parser.parse_args()

import sys
import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import networkx as nx
import errno    
import os
import math
from collections import Counter
from Bio.Alphabet import generic_dna

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def logger(ltype,msg):
    st = datetime.datetime.now()
    if (ltype=="err" or ltype=="warn"):
        sys.stderr.write(" :: ".join([str(st),ltype,msg])+"\n")
        if ltype=="err":
            sys.exit(1)
    elif ltype=="info" or ltype=="debug":
        sys.stdout.write(" :: ".join([str(st),ltype,msg])+"\n")
    else:
        logger("warn","error in the code the logger function is called with a base type...")

rna_gap=int(args.rna_gap)

if args.v:
    logger("info","Begin the analysis with Scaff2link version "+version)

mkdir_p(args.outDir)

# pour avoir les arguments c'est args.fasta / args.phylo / args.rna etc ... ça te renvoie la string (type string)

g=nx.MultiDiGraph()

## Adding all vertexes to the graph
seqDict=dict()

for record in SeqIO.parse(args.fasta, "fasta"):
    seqDict[record.id]=record.seq
    g.add_node(record.id,
    length=len(record.seq),
    log10Len=math.log10(len(record.seq)),
    strand='')

if args.v : logger("info","Number of nodes parsed: "+str(len(g.nodes)))

# Adding edges from Ragout :

if args.v : logger("info","Parsing edges from synteny information using ragout: "+args.phylo)

scaffolding_type_1="synteny"
current_scaffold=""
fromName=""
gap_length=0
with open(args.phylo) as f:
    for line in f:
        line=line.rstrip("\n")
        if line[0] != "#" and line !="":
            elems=line.split("\t")
            if current_scaffold == elems[0]:
                if elems[4] == "N":
                    gap_length=int(elems[5])
                else:
                    if not (fromName in g.nodes and elems[5] in g.nodes):
                        logger("err","Unknown node name when parsing ragout edge from "+fromName+" to "+elems[5])
                    g.add_edge(fromName,elems[5],type=scaffolding_type_1,fromStrand=from_orientation,toStrand=elems[8],gap=gap_length, readsCount=-1)
                    g.nodes[fromName]["strand"]=from_orientation
                    g.nodes[elems[5]]["strand"]=elems[8]
                    fromName=elems[5]
                    from_orientation=elems[8]
            else:
                current_scaffold=elems[0]
                fromName=elems[5]
                from_orientation=elems[8]
f.close()

if args.v : logger("info","Parsing joint pairs in edges from rna-seq information using agouti: "+args.rna)

def strand_reverse(strand):
    if strand=="+":
        sr="-"
    elif strand=="-":
        sr="+"
    else:
        logger("err","One strand information is not + or -, please check input files")
    return(sr)

scaffolding_type_2="rna-seq"

joint_pairs=dict()

with open(args.rna) as f:
    for line in f:
        line=line.rstrip()
        elems=line.split("\t")
        fromName=elems[1]
        toName=elems[4]
        if args.lib_type[0]=='f':
            FromStrand=elems[3]
        else:
            FromStrand=strand_reverse(elems[3])
        if args.lib_type[1]=='f':
            ToStrand=elems[6]
        else:
            ToStrand=strand_reverse(elems[6])



        key="\t".join([FromStrand,fromName,ToStrand,toName])
        if key not in joint_pairs.keys():
            joint_pairs[key]=1
        else:
            joint_pairs[key]+=1

def add_stranded_edge(fromName,toName,FromStrand,ToStrand,scaffolding_type,readsCount):
    global g
    edgeToAdd=False
    consistencyIssue=False
    if (g.nodes[fromName]["strand"]==""):
        if g.nodes[toName]["strand"]=="":
            g.add_edge(fromName,toName,type=scaffolding_type,fromStrand=FromStrand,toStrand=ToStrand,gap=-1,readsCount=readsCount)
            g.nodes[fromName]["strand"]=FromStrand
            g.nodes[toName]["strand"]=ToStrand
        elif g.nodes[toName]["strand"]==".":
            g.add_edge(fromName,toName,type=scaffolding_type,fromStrand=FromStrand,toStrand=ToStrand,gap=-1,readsCount=readsCount)
            g.nodes[fromName]["strand"]=FromStrand
        elif g.nodes[toName]["strand"]==ToStrand:
            g.add_edge(fromName,toName,type=scaffolding_type,fromStrand=FromStrand,toStrand=ToStrand,gap=-1,readsCount=readsCount)
            g.nodes[fromName]["strand"]=FromStrand
        elif g.nodes[toName]["strand"]==strand_reverse(ToStrand):
            g.add_edge(toName,fromName,type=scaffolding_type,fromStrand=strand_reverse(ToStrand),toStrand=strand_reverse(FromStrand),gap=-1,readsCount=readsCount)
            g.nodes[fromName]["strand"]=strand_reverse(FromStrand)
        else:
            logger("err","A case was not anticipated to add stranded edge with undefined strand on the coming node")
    elif (g.nodes[toName]["strand"]==""):
        if g.nodes[fromName]["strand"]==".":
            g.add_edge(fromName,toName,type=scaffolding_type,fromStrand=FromStrand,toStrand=ToStrand,gap=-1,readsCount=readsCount)
            g.nodes[toName]["strand"]=ToStrand
        elif g.nodes[fromName]["strand"]==FromStrand:
            g.add_edge(fromName,toName,type=scaffolding_type,fromStrand=FromStrand,toStrand=ToStrand,gap=-1,readsCount=readsCount)
            g.nodes[toName]["strand"]=ToStrand
        elif g.nodes[fromName]["strand"]==strand_reverse(FromStrand):
            g.add_edge(toName,fromName,type=scaffolding_type,fromStrand=strand_reverse(ToStrand),toStrand=strand_reverse(FromStrand),gap=-1,readsCount=readsCount)
            g.nodes[toName]["strand"]=strand_reverse(ToStrand)
        else:
            logger("err","A case was not anticipated to add stranded edge with undefined strand on the arrival node")
    elif (g.nodes[fromName]["strand"]==FromStrand) and (g.nodes[toName]["strand"]==ToStrand):
        numFromTo=g.number_of_edges(fromName,toName)
        numToFrom=g.number_of_edges(toName,fromName)
        edgeToAdd=True
        reverseAdd=False
    elif (g.nodes[fromName]["strand"]==strand_reverse(FromStrand)) and (g.nodes[toName]["strand"]==strand_reverse(ToStrand)):
        numFromTo=g.number_of_edges(fromName,toName)
        numToFrom=g.number_of_edges(toName,fromName)
        edgeToAdd=True
        reverseAdd=True

    if (g.nodes[fromName]["strand"]=="."):
        numFromTo=g.number_of_edges(fromName,toName)
        numToFrom=g.number_of_edges(toName,fromName)
        if g.nodes[toName]["strand"]==".":
            if (numFromTo+numToFrom)==0:
                g.add_edge(fromName,toName,type=scaffolding_type,fromStrand=FromStrand,toStrand=ToStrand,gap=-1,readsCount=readsCount)
            else:
                edgeToAdd=True
                reverseAdd=False
        elif g.nodes[toName]["strand"]==ToStrand:
            if (numFromTo+numToFrom)==0:
                g.add_edge(fromName,toName,type=scaffolding_type,fromStrand=FromStrand,toStrand=ToStrand,gap=-1,readsCount=readsCount)
            else:
                edgeToAdd=True
                reverseAdd=False
        elif g.nodes[toName]["strand"]==strand_reverse(ToStrand):
            if (numFromTo+numToFrom)==0:
                g.add_edge(toName,fromName,type=scaffolding_type,fromStrand=strand_reverse(ToStrand),toStrand=strand_reverse(FromStrand),gap=-1,readsCount=readsCount)
            else:
                edgeToAdd=True
                reverseAdd=True
        else:
            logger("err","A case was not anticipated to add stranded edge with unconsitant strand on coming node")
    elif (g.nodes[toName]["strand"]=="."):
        numFromTo=g.number_of_edges(fromName,toName)
        numToFrom=g.number_of_edges(toName,fromName)
        if g.nodes[fromName]["strand"]==FromStrand:
            if (numFromTo+numToFrom)==0:
                g.add_edge(fromName,toName,type=scaffolding_type,fromStrand=FromStrand,toStrand=ToStrand,gap=-1,readsCount=readsCount)
            else:
                edgeToAdd=True
                reverseAdd=False
        elif g.nodes[fromName]["strand"]==strand_reverse(FromStrand):
            if (numFromTo+numToFrom)==0:
                g.add_edge(toName,fromName,type=scaffolding_type,fromStrand=strand_reverse(ToStrand),toStrand=strand_reverse(FromStrand),gap=-1,readsCount=readsCount)
            else:
                edgeToAdd=True
                reverseAdd=True
        else:
            logger("err","A case was not anticipated to add stranded edge with unconsitant strand on arrival node")

    if edgeToAdd:
            if numFromTo>0:
                for i in range(0,numFromTo):
                    if ((g.edges[fromName,toName,i]['fromStrand']==FromStrand) and (g.edges[fromName,toName,i]['toStrand']==ToStrand)):
                        g.edges[fromName,toName,i]["type"]=g.edges[fromName,toName,i]["type"]+'_AND_'+scaffolding_type
                        g.edges[fromName,toName,i]['readsCount']=readsCount
                        edgeToAdd=False
            if numToFrom>0:
                for i in range(0,numToFrom):
                    # Orientation consistency ?
                    if ((g.edges[toName,fromName,i]['fromStrand']==strand_reverse(ToStrand)) and (g.edges[toName,fromName,i]['toStrand']==strand_reverse(FromStrand))):
                        g.edges[toName,fromName,i]["type"]=scaffolding_type_1+'_AND_'+scaffolding_type_2
                        g.edges[toName,fromName,i]['readsCount']=joint_pairs[key]
                        edgeToAdd=False
            if edgeToAdd:
                if reverseAdd:
                    g.add_edge(toName,fromName,type=scaffolding_type,fromStrand=strand_reverse(ToStrand),toStrand=strand_reverse(FromStrand),gap=-1,readsCount=readsCount)
                else:
                    g.add_edge(fromName,toName,type=scaffolding_type,fromStrand=FromStrand,toStrand=ToStrand,gap=-1,readsCount=readsCount)

if args.v : logger("info","Adding rna-seq based edges to the scaffolding graph")

for key in joint_pairs.keys():
    if joint_pairs[key]>(int(args.min_reads)-1):
        key_list=key.split("\t")
        FromStrand=key_list[0]
        fromName=key_list[1]
        ToStrand=key_list[2]
        toName=key_list[3]
        if not (fromName in g.nodes and toName in g.nodes):
            logger("err","Unknown node name when parsing agouti edge from "+fromName+" to "+toName)
        # Here we add the edge, but try to simplify it if it already exist because of Ragout
        numFromTo=g.number_of_edges(fromName,toName)
        add_stranded_edge(fromName,toName,FromStrand,ToStrand,scaffolding_type_2,joint_pairs[key])

if args.v : logger("info","Scaffolding graph completed, collecting nodes statistics")

TotalNodes=len(g.nodes)
ConnectedNodes=TotalNodes
ConsistantNodes=TotalNodes

node2remove=list()

for node in g.nodes:
    if nx.is_isolate(g,node): # WAY TO GET THE DEGREE OF NODE
        #g.remove_node(node)
        node2remove.append(node)
        ConnectedNodes-=1
        ConsistantNodes-=1
    elif g.nodes[node]["strand"]=='.':
        ConsistantNodes-=1

for node in node2remove:
    g.remove_node(node)

logger("info","Out of "+str(TotalNodes)+" initial nodes, "+str(ConnectedNodes)+" are connected, and among them "+str(ConsistantNodes)+" are strand consistent")

if args.v : 
    logger("info","Collecting edges statistics")
    edgeCounter=Counter(list(g.edges))
    t1=0
    t2=0
    t12=0
    edges_list=list(g.edges)
    edges_set=set(g.edges)
    for edge in edgeCounter.keys():
        for key in range(0,edgeCounter[edge]):
            if g[edge[0]][edge[1]][key]["type"]==scaffolding_type_1:
                t1+=1
            elif g[edge[0]][edge[1]][key]["type"]==scaffolding_type_2:
                t2+=1
            elif g[edge[0]][edge[1]][key]["type"]==scaffolding_type_1+"_AND_"+scaffolding_type_2:
                t12+=1
            else:
                logger("err","unknown scaffold type")

logger("info","Edges statistics : \n"+
    ";".join([scaffolding_type_1,scaffolding_type_2,scaffolding_type_1+"_AND_"+scaffolding_type_2])+
    "\n"+";".join([str(t1),str(t2),str(t12)]))

mkdir_p(args.outDir+'/00-complete_graph')

nx.write_graphml(g, args.outDir+"/00-complete_graph/graph.graphml")

if args.v : logger("info","starting first chain simplification")

def chain_simplification(chain,start,end,name):
    global g
    global seqDict
    global rna_gap
    def_chain=start
    chain_g=g.subgraph(chain)
    setchain=set(chain)
    setchain.remove(start)

    if chain_g.nodes[start]["strand"]=="+":
        seqDict[name]=seqDict[start]
    elif chain_g.nodes[start]["strand"]=="-":
        seqDict[name]=seqDict[start].reverse_complement()
    else:
        logger("err","unexpected strand value at initiation of chain simplification")
    del seqDict[start]
    node=start

    while len(setchain)!=0:
        lastnode=node
        node=list(chain_g.out_edges(node))[0][1]
        setchain.remove(node)
        if chain_g[lastnode][node][0]["gap"]==-1:
            gap=rna_gap*'n'
            gap_value=rna_gap
            def_chain+='('+str(gap_value)+"n)"+chain_g.nodes[node]["strand"]+node
        else:
            gap=chain_g[lastnode][node][0]["gap"]*'N'
            gap_value=chain_g[lastnode][node][0]["gap"]
            def_chain+='('+str(gap_value)+"N)"+chain_g.nodes[node]["strand"]+node
        if chain_g.nodes[node]["strand"]=="+":
            seqDict[name]+=Seq(gap,generic_dna)+seqDict[node]
        elif chain_g.nodes[node]["strand"]=="-":
            seqDict[name]+=Seq(gap,generic_dna)+seqDict[node].reverse_complement()
        else:
            logger("err","unexpected strand value during chain simplification")
        del seqDict[node]
    inEdges=list(g.in_edges(start))
    outEdges=list(g.out_edges(end))

    g.add_node(name,
        length=len(seqDict[name]),
        log10Len=math.log10(len(seqDict[name])),
        strand='+')

    for inEdge in set(inEdges):
        for key in range(0,inEdges.count(inEdge)):
            g.add_edge(inEdge[0],name,
                type=g[inEdge[0]][inEdge[1]][key]["type"],
                fromStrand=g[inEdge[0]][inEdge[1]][key]["fromStrand"],
                toStrand='+',
                gap=g[inEdge[0]][inEdge[1]][key]["gap"],
                readsCount=g[inEdge[0]][inEdge[1]][key]["readsCount"]
                )
    for outEdge in set(outEdges):
        for key in range(0,outEdges.count(outEdge)):
            g.add_edge(name,outEdge[1],
                type=g[outEdge[0]][outEdge[1]][key]["type"],
                fromStrand="+",
                toStrand=g[outEdge[0]][outEdge[1]][key]["toStrand"],
                gap=g[outEdge[0]][outEdge[1]][key]["gap"],
                readsCount=g[outEdge[0]][outEdge[1]][key]["readsCount"]
                )
    g.remove_nodes_from(chain)

    return(def_chain)

setNodes=set(g.nodes)

basename="scaff2links_chain_"
count_name=1

while len(setNodes)!=0:
    try:
        node=setNodes.pop()
        if g.nodes[node]["strand"]!='.':
            inspectIn=True
            inspectOut=True
            chain=[node]
            nodeIn=node
            nodeOut=node
            if args.v: logger('info',"Chain simplification started on node:"+node)
            while inspectIn:
                if len(g.in_edges(nodeOut))==1 :
                    possibleNodeOut=list(g.in_edges(nodeOut))[0][0]
                    if args.d: logger('debug',"IN chain extension: (in)"+nodeOut+"\t(out)"+possibleNodeOut)
                    if (len(g.out_edges(possibleNodeOut))==1) and (g.nodes[possibleNodeOut]["strand"]!='.'):
                        if args.d: logger('debug',"CONTINUE")
                        nodeOut=possibleNodeOut
                        chain.append(nodeOut)
                        setNodes.remove(nodeOut)
                    else : 
                        inspectIn=False
                        if args.d: logger('debug',"STOP")
                else: inspectIn=False
            while inspectOut:
                if len(g.out_edges(nodeIn))==1 :
                    possibleNodeIn=list(g.out_edges(nodeIn))[0][1]
                    if args.d: logger('debug',"OUT chain extension: (out)"+nodeIn+"\t(in)"+possibleNodeIn)
                    if (len(g.in_edges(possibleNodeIn))==1) and (g.nodes[possibleNodeIn]["strand"]!='.'):
                        if args.d: logger('debug',"CONTINUE")
                        nodeIn=possibleNodeIn
                        chain.append(nodeIn)
                        setNodes.remove(nodeIn)
                    else: 
                        inspectOut=False
                        if args.d: logger('debug',"STOP")
                else: inspectOut=False
            name=basename+str(count_name)
            count_name+=1
            def_chain=chain_simplification(chain=chain,start=nodeOut,end=nodeIn,name=name)
            if args.d: logger('debug',name+"\t"+def_chain)
    except:
        nx.write_graphml(g, args.outDir+"/error.graphml")
        logger("err","error during first chain simplification, the current graph have been written")

mkdir_p(args.outDir+'/01-simplified_graph')

nx.write_graphml(g, args.outDir+"/01-simplified_graph/graph.graphml")


with open(args.outDir+"/01-simplified_graph/scaffolds.fasta", "w") as output_handle:
    for key in seqDict.keys():
        SeqIO.write(SeqRecord(seqDict[key],id=key,description=''), output_handle, "fasta")