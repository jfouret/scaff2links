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
parser.add_argument('-log', metavar='/path', required=False,default='./scaff2links.log', help="""log file""")
parser.add_argument('-outDir', metavar='/path', required=True, help="""Output directory""")
parser.add_argument('-v',action='store_true', help="verbose")
parser.add_argument('-d',action='store_true', help="Debug info")
args=parser.parse_args()

import sys
import datetime
import traceback
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import networkx as nx
import errno    
import os
import math
import logging
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
def strand_reverse(strand):
    if strand=="+":
        sr="-"
    elif strand=="-":
        sr="+"
    else:
        logging.critical("One strand information is not + or -, please check input files")
        sys.exit(1)
    return(sr)
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
            logging.critical("A case was not anticipated to add stranded edge with undefined strand on the coming node")
            sys.exit(1)
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
            logging.critical("A case was not anticipated to add stranded edge with undefined strand on the arrival node")
            sys.exit(1)
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
            logging.critical("A case was not anticipated to add stranded edge with unconsitant strand on coming node")
            sys.exit(1)
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
            logging.critical("A case was not anticipated to add stranded edge with unconsitant strand on arrival node")
            sys.exit(1)

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
                        g.edges[toName,fromName,i]['readsCount']=readsCount
                        edgeToAdd=False
            if edgeToAdd:
                if reverseAdd:
                    g.add_edge(toName,fromName,type=scaffolding_type,fromStrand=strand_reverse(ToStrand),toStrand=strand_reverse(FromStrand),gap=-1,readsCount=readsCount)
                else:
                    g.add_edge(fromName,toName,type=scaffolding_type,fromStrand=FromStrand,toStrand=ToStrand,gap=-1,readsCount=readsCount)
def chain_simplification(chain,start,end,name,chain_is_path):
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
        logging.critical("unexpected strand value at initiation of chain simplification")
        sys.exit(1)
    del seqDict[start]
    node=start
    while len(setchain)!=0:
        lastnode=node
        if chain_is_path:
            node=chain[len(chain)-len(setchain)]
        else:
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
            logging.critical("unexpected strand value during chain simplification")
            sys.exit(1)
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

def main():
    global seqDict
    global rna_gap
    global g
    global args
    global scaffolding_type_1
    global scaffolding_type_2

    """
    Parsing arguments and initiating variables
    """

    if args.d:
        logging.basicConfig(format='%(asctime)s::%(levelname)s::%(message)s',filename=args.log, level=logging.DEBUG)
    elif args.v:
        logging.basicConfig(format='%(asctime)s::%(levelname)s::%(message)s',filename=args.log, level=logging.INFO)
    else:
        logging.basicConfig(format='%(asctime)s::%(levelname)s::%(message)s',filename=args.log)
    rna_gap=int(args.rna_gap)
    logging.info("Begin the analysis with Scaff2link version "+version)
    mkdir_p(args.outDir)
    # pour avoir les arguments c'est args.fasta / args.phylo / args.rna etc ... ça te renvoie la string (type string)
    g=nx.MultiDiGraph()
    ## Adding all vertexes to the graph

    """
    Reading fasta
    """

    seqDict=dict()
    for record in SeqIO.parse(args.fasta, "fasta"):
        seqDict[record.id]=record.seq
        g.add_node(record.id,
        length=len(record.seq),
        log10Len=math.log10(len(record.seq)),
        strand='')
    logging.info("Number of nodes parsed: "+str(len(g.nodes)))

    """
    Adding edges from Ragout :
    """

    logging.info("Parsing edges from synteny information using ragout: "+args.phylo)
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
                            logging.critical("Unknown node name when parsing ragout edge from "+fromName+" to "+elems[5])
                            sys.exit(1)
                        g.add_edge(fromName,elems[5],type=scaffolding_type_1,fromStrand=from_orientation,toStrand=elems[8],gap=gap_length, readsCount=-1)
                        g.nodes[fromName]["strand"]=from_orientation
                        g.nodes[elems[5]]["strand"]=elems[8]
                        fromName=elems[5]
                        from_orientation=elems[8]
                else:
                    current_scaffold=elems[0]
                    fromName=elems[5]
                    from_orientation=elems[8]

    """
    Parsing edges from Agouti :
    """

    logging.info("Parsing joint pairs in edges from rna-seq information using agouti: "+args.rna)
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

    """
    Adding parsed edges from Agouti :
    """

    logging.info("Adding rna-seq based edges to the scaffolding graph")
    for key in joint_pairs.keys():
        if joint_pairs[key]>(int(args.min_reads)-1):
            key_list=key.split("\t")
            FromStrand=key_list[0]
            fromName=key_list[1]
            ToStrand=key_list[2]
            toName=key_list[3]
            if not (fromName in g.nodes and toName in g.nodes):
                logging.critical("Unknown node name when parsing agouti edge from "+fromName+" to "+toName)
                sys.exit(1)
            # Here we add the edge, but try to simplify it if it already exist because of Ragout
            numFromTo=g.number_of_edges(fromName,toName)
            add_stranded_edge(fromName,toName,FromStrand,ToStrand,scaffolding_type_2,joint_pairs[key])

    """
    Finalizing and reporting the initial graph (end of step 00)
    """

    logging.info("Scaffolding graph completed, collecting nodes statistics")
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
    logging.info("Out of "+str(TotalNodes)+" initial nodes, "+str(ConnectedNodes)+" are connected, and among them "+str(ConsistantNodes)+" are strand consistent")
    if args.v : 
        logging.info("Collecting edges statistics")
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
                    logging.critical("unknown scaffold type")
                    sys.exit(1)
        logging.info("Edges statistics : \n"+
            ";".join([scaffolding_type_1,scaffolding_type_2,scaffolding_type_1+"_AND_"+scaffolding_type_2])+
            "\n"+";".join([str(t1),str(t2),str(t12)]))
    mkdir_p(args.outDir+'/00-complete_graph')
    nx.write_graphml(g, args.outDir+"/00-complete_graph/graph.graphml")

    """
    Step 01 : chain simplification
    """

    logging.info("starting first chain simplification")
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
                logging.debug("Chain simplification started on node:"+node)
                while inspectIn:
                    if len(g.in_edges(nodeOut))==1 :
                        possibleNodeOut=list(g.in_edges(nodeOut))[0][0]
                        if g.has_edge(nodeOut,possibleNodeOut):
                            inspectIn=False
                            logging.debug("STOP")
                        else:
                            logging.debug("IN chain extension: (in)"+nodeOut+"\t(out)"+possibleNodeOut)
                            if (len(g.out_edges(possibleNodeOut))==1) and (g.nodes[possibleNodeOut]["strand"]!='.'):
                                logging.debug("CONTINUE")
                                nodeOut=possibleNodeOut
                                chain.append(nodeOut)
                                setNodes.remove(nodeOut)
                            else : 
                                inspectIn=False
                                logging.debug("STOP")
                    else: inspectIn=False
                while inspectOut:
                    if len(g.out_edges(nodeIn))==1 :
                        possibleNodeIn=list(g.out_edges(nodeIn))[0][1]
                        if g.has_edge(possibleNodeIn,nodeIn):
                            inspectOut=False
                            logging.debug("STOP")
                        else:
                            logging.debug("OUT chain extension: (out)"+nodeIn+"\t(in)"+possibleNodeIn)
                            if (len(g.in_edges(possibleNodeIn))==1) and (g.nodes[possibleNodeIn]["strand"]!='.'):
                                logging.debug("CONTINUE")
                                nodeIn=possibleNodeIn
                                chain.append(nodeIn)
                                setNodes.remove(nodeIn)
                            else: 
                                inspectOut=False
                                logging.debug("STOP")
                    else: inspectOut=False
                name=basename+str(count_name)
                count_name+=1
                def_chain=chain_simplification(chain=chain,start=nodeOut,end=nodeIn,name=name,chain_is_path=False)
                logging.debug(name+"\t"+def_chain)
        except:
            nx.write_graphml(g, args.outDir+"/error.graphml")
            print(traceback.format_exc())
            logging.critical("error during first chain simplification, the current graph have been written")
            sys.exit(1)

    """
    end of step 01 writing...
    """

    mkdir_p(args.outDir+'/01-simplified_graph')
    nx.write_graphml(g, args.outDir+"/01-simplified_graph/graph.graphml")
    node2remove=list()
    for node in g.nodes:
        if nx.is_isolate(g,node): # WAY TO GET THE DEGREE OF NODE
            #g.remove_node(node)
            node2remove.append(node)
    for node in node2remove:
        g.remove_node(node)
    with open(args.outDir+"/01-simplified_graph/scaffolds.fasta", "w") as output_handle:
        for key in seqDict.keys():
            SeqIO.write(SeqRecord(seqDict[key],id=key,description=''), output_handle, "fasta")

    """
    Step 02 : dag simplification
    """

    logging.info("starting dag simplification")
    logging.info("Find remaining bridges")
    g_broken=nx.Graph(g.copy())
    edge2remove=list()
    for e in nx.bridges(g_broken):
        edge2remove.append(e)
    for e_ in edge2remove:
        g_broken.remove_edge(e_[0],e_[1])
    basename="scaff2links_dag_"
    count_name=1
    logging.info("Connected component analysis for DAG")
    for nset in nx.connected_components(g_broken):
        subg=g.subgraph(nset)
        if nx.is_directed_acyclic_graph(subg) and len(list(subg.nodes()))>1:
            lpath=nx.dag_longest_path(subg)
            logging.debug("DAG found :"+str(nset)+"| type="+str(type(nset)))
            if len(list(subg.nodes()))==2:
                logging.warn("A connected component resulting from the graph where all bridges have been removed is of size 2, which is mathematically unexpected")
            if len(lpath)==len(nset):
                logging.debug("DAG with all nodes in the longest path: ("+")->-(".join(lpath)+')')
                # check neighbor
                pos=-1
                start=0
                for node in lpath:
                    pos+=1
                    add_set=set()
                    if node==lpath[0] or pos == start:
                        for e in list(g.in_edges(node)):
                            add_set.add(e[0])
                    if node==lpath[-1]:
                        for e in list(g.out_edges(node)):
                            add_set.add(e[1])
                    if not (set(nx.all_neighbors(g,node)) <= (add_set | nset)):
                        logging.debug("Found cutting node in DAG: "+lpath[pos])
                        for e in list(g.out_edges(node)):
                            add_set.add(e[1])
                        if (set(nx.all_neighbors(g,node)) <= (add_set | nset)) and (pos-start>0):
                            logging.debug("cutted DAG simplification (including the cutting node) at "+lpath[pos]+": ("+")->-(".join([ lpath[x] for x in range(start,pos+1) ] )+')')
                            chain_simplification([ lpath[x] for x in range(start,pos+1) ],lpath[start],lpath[pos],basename+str(count_name),chain_is_path=True)
                            count_name+=1
                            start=pos+1
                            checkStart=False
                        elif pos-start>1:
                            logging.debug("cutted DAG simplification (excluding the cutting node) at "+lpath[pos]+": ("+")->-(".join([ lpath[x] for x in range(start,pos) ] )+')')
                            chain_simplification([ lpath[x] for x in range(start,pos) ],lpath[start],lpath[pos-1],basename+str(count_name),chain_is_path=True)
                            count_name+=1
                            checkStart=True
                        else: checkStart=True
                        if checkStart:
                            start=pos
                            for e in list(g.out_edges(node)):
                                if not e[1] in nset:
                                    start=pos+1
                                    logging.debug("NB: Cutting node excluded as start")
                                    break
                            if start==pos:
                                logging.debug("NB: Cutting node included as start")
                if pos-start>0:
                    logging.debug("Final DAG simplification : ("+")->-(".join([ lpath[x] for x in range(start,pos+1) ] )+')')
                    chain_simplification([ lpath[x] for x in range(start,pos+1) ],lpath[start],lpath[pos],basename+str(count_name),chain_is_path=True)
                    count_name+=1
                start=pos+1

    """
    end of step 02 writing...
    """

    logging.info("Write graph and fasta after DAG simplification")
    mkdir_p(args.outDir+'/02-after_dag_graph')
    nx.write_graphml(g, args.outDir+"/02-after_dag_graph/graph.graphml")
    node2remove=list()
    for node in g.nodes:
        if nx.is_isolate(g,node): # WAY TO GET THE DEGREE OF NODE
            #g.remove_node(node)
            node2remove.append(node)
    for node in node2remove:
        g.remove_node(node)
    with open(args.outDir+"/02-after_dag_graph/scaffolds.fasta", "w") as output_handle:
        for key in seqDict.keys():
            SeqIO.write(SeqRecord(seqDict[key],id=key,description=''), output_handle, "fasta")

try:
    main()
except:
    logging.critical(traceback.format_exc())
    sys.exit(1)
sys.exit(0)