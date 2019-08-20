import os
import random
import json
import geopandas as gpd
import functools
import datetime
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv
from networkx.readwrite import json_graph
import math
import seaborn as sns
from functools import partial
import networkx as nx
import numpy as np


from sympy.matrices import Matrix


from gerrychain import Graph
from gerrychain import MarkovChain
from gerrychain.constraints import (Validator, single_flip_contiguous,
within_percent_of_ideal_population, UpperBound)
from gerrychain.proposals import propose_random_flip, propose_chunk_flip
from gerrychain.accept import always_accept
from gerrychain.updaters import Election,Tally,cut_edges
from gerrychain import GeographicPartition
from gerrychain.partition import Partition
from gerrychain.proposals import recom
from gerrychain.metrics import mean_median, efficiency_gap


mu = 2.63815853
bases = [1]#[.1,1/mu**2,.2,1/mu,.8,1,mu,4,mu**2,10]
base=.1
pops=[.05]


def fixed_endpoints(partition):
    return partition.assignment[(19,0)] != partition.assignment[(20,0)] and partition.assignment[(19,39)] != partition.assignment[(20,39)]
    
    
def boundary_condition(partition):
    
    blist = partition["boundary"]
    o_part = partition.assignment[blist[0]]

    for x in blist:
        if partition.assignment[x] != o_part:
            return True
    
    return False
    
    
def boundary_slope(partition):

    a=[]
    b=[]
    c=[]
    d=[]
    e=[]
    
    for x in partition["cut_edges"]:
        if x[0][0] == 0 and x[1][0] == 0:
            a.append(x)
        elif x[0][1] == 0 and x[1][1] == 0:
            b.append(x)
        elif x[0][0] == 39 and x[1][0] == 39:
            c.append(x)
        elif x[0][1] == 39 and x[1][1] == 39:
            d.append(x)
        elif x in [((0,1),(1,0)), ((0,38),(1,39)), ((38,0),(39,1)), ((38,39),(39,38))]:
            e.append(x)
        elif x in [((1,0),(0,1)), ((1,39),(0,38)), ((39,1),(38,0)), ((39,38),(38,39))]:
            e.append(x)
                    

    return list(set(a+b+c+d+e))
    

def annealing_cut_accept_backwards(partition):
    boundaries1  = {x[0] for x in partition["cut_edges"]}.union({x[1] for x in partition["cut_edges"]})
    boundaries2  = {x[0] for x in partition.parent["cut_edges"]}.union({x[1] for x in partition.parent["cut_edges"]})
    
    t = partition["step_num"]
    
    
    #if t <100000:
    #    beta = 0  
    #elif t<400000:
    #    beta = (t-100000)/100000 #was 50000)/50000
    #else:
    #    beta = 3
    base = .1
    beta = 5
        
    bound = 1
    if partition.parent is not None:
        bound = (base**(beta*(-len(partition["cut_edges"])+len(partition.parent["cut_edges"]))))*(len(boundaries1)/len(boundaries2))  
        
        if not popbound(partition):
            bound = 0         
        if not single_flip_contiguous(partition):
            bound = 0             
        #bound = min(1, (how_many_seats_value(partition, col1="G17RATG",
         #col2="G17DATG")/how_many_seats_value(partition.parent, col1="G17RATG",
         #col2="G17DATG"))**2  ) #for some states/elections probably want to add 1 to denominator so you don't divide by zero


    return random.random() < bound
    
    
def go_nowhere(partition):
    return partition.flip(dict())
    
    
def slow_reversible_propose(partition):
    """Proposes a random boundary flip from the partition in a reversible fasion
    by selecting uniformly from the (node, flip) pairs.
    Temporary version until we make an updater for this set.
    :param partition: The current partition to propose a flip from.
    :return: a proposed next `~gerrychain.Partition`
    """

    #b_nodes = {(x[0], partition.assignment[x[1]]) for x in partition["cut_edges"]
    #           }.union({(x[1], partition.assignment[x[0]]) for x in partition["cut_edges"]})

    flip = random.choice(list(partition["b_nodes"]))
    
    return partition.flip({flip[0]: flip[1]})

def slow_reversible_propose_bi(partition):
    """Proposes a random boundary flip from the partition in a reversible fasion
    by selecting uniformly from the (node, flip) pairs.
    Temporary version until we make an updater for this set.
    :param partition: The current partition to propose a flip from.
    :return: a proposed next `~gerrychain.Partition`
    """

    #b_nodes = {(x[0], partition.assignment[x[1]]) for x in partition["cut_edges"]
    #           }.union({(x[1], partition.assignment[x[0]]) for x in partition["cut_edges"]})

    fnode = random.choice(list(partition["b_nodes"]))
    
    return partition.flip({fnode: -1*partition.assignment[fnode]})

def geom_wait(partition):
    return int(np.random.geometric(len(list(partition["b_nodes"]))/(len(partition.graph.nodes)**(len(partition.parts))-1) ,1))-1
    

def b_nodes(partition):
    return {(x[0], partition.assignment[x[1]]) for x in partition["cut_edges"]
               }.union({(x[1], partition.assignment[x[0]]) for x in partition["cut_edges"]})

def b_nodes_bi(partition):
    return {x[0] for x in partition["cut_edges"]}.union({x[1] for x in partition["cut_edges"]})   


def uniform_accept(partition):

    bound = 0
    if popbound(partition) and single_flip_contiguous(partition) and boundary_condition(partition):
        bound = 1
        
    return random.random() < bound
        
#BUILD GRAPH



def cut_accept(partition):
      
    bound = 1
    if partition.parent is not None:
        bound = (partition["base"]**(-len(partition["cut_edges"])+len(partition.parent["cut_edges"])))#*(len(boundaries1)/len(boundaries2))
 
       

    return random.random() < bound


for pop1 in pops:
    for base in bases:
        for alignment in [0]:

            gn=25
            k=2
            ns=100
            p=.5
    
            graph=nx.grid_graph([k*gn,k*gn])
    
            ########## BUILD ASSIGNMENT
            #cddict = {x: int(x[0]/gn)  for x in graph.nodes()}        
            cddict = {x: 1-2*int(x[0]/gn)  for x in graph.nodes()}
            
    
    
            for n in graph.nodes():
                graph.node[n]["population"]=1
                graph.node[n]["part_sum"]=cddict[n]
                graph.node[n]["last_flipped"]=0
                graph.node[n]["num_flips"]=0
               
                if random.random()<p:
                    graph.node[n]["pink"]=1
                    graph.node[n]["purple"]=0
                else:
                    graph.node[n]["pink"]=0
                    graph.node[n]["purple"]=1
                if 0 in n or k*gn-1 in n:
                    graph.node[n]["boundary_node"]=True
                    graph.node[n]["boundary_perim"]=1
    
                else:
                    graph.node[n]["boundary_node"]=False


            for edge in graph.edges():
                graph[edge[0]][edge[1]]['cut_times'] = 0 
    
            #this part adds queen adjacency
            #for i in range(k*gn-1):
            #    for j in range(k*gn):
            #        if j<(k*gn-1):
            #            graph.add_edge((i,j),(i+1,j+1))
            #            graph[(i,j)][(i+1,j+1)]["shared_perim"]=0
            #        if j >0:
            #            graph.add_edge((i,j),(i+1,j-1))
            #            graph[(i,j)][(i+1,j-1)]["shared_perim"]=0
                    
                    
            ######PLOT GRIDS
            """      
            plt.figure()
            nx.draw(graph, pos = {x:x for x in graph.nodes()} ,node_size = ns, node_shape ='s')
            plt.show()
    
            cdict = {1:'pink',0:'purple'}
    
            plt.figure()
            nx.draw(graph, pos = {x:x for x in graph.nodes()}, node_color = [cdict[graph.node[x]["pink"]] for x in graph.nodes()],node_size = ns, node_shape ='s' )
            plt.show()
    
            plt.figure()
            nx.draw(graph, pos = {x:x for x in graph.nodes()}, node_color = [cddict[x] for x in graph.nodes()] ,node_size = ns, node_shape ='s',cmap = 'tab20')
            plt.show()
            """
            ####CONFIGURE UPDATERS
            

            def step_num(partition):
                parent = partition.parent
    
                if not parent:
                    return 0
    
    
                return parent["step_num"] + 1
                
                
            bnodes = [x for x in graph.nodes() if graph.node[x]["boundary_node"] ==1]
    
            def bnodes_p(partition):
    
                
                return [x for x in graph.nodes() if graph.node[x]["boundary_node"] ==1]
    
            updaters = {'population': Tally('population'),
                                #"boundary":bnodes_p,
                                "slope": boundary_slope,
                                'cut_edges': cut_edges,
                                'step_num': step_num,
                                'b_nodes':b_nodes_bi,
                                #'base':new_base,
                                #'geom':geom_wait,
                                #"Pink-Purple": Election("Pink-Purple", {"Pink":"pink","Purple":"purple"})
                    }                  
                                
                                
    
                                
    
            #########BUILD PARTITION
    
            grid_partition = Partition(graph,assignment=cddict,updaters=updaters)
            
            
            colors = {1:'darkorange', -1:'mediumseagreen'}
    
            #ADD CONSTRAINTS
            popbound=within_percent_of_ideal_population(grid_partition,pop1) 

            plt.figure()
            nx.draw(graph, pos = {x:x for x in graph.nodes()}, node_color = [colors[dict(grid_partition.assignment)[x]] for x in graph.nodes()] ,node_size = ns, node_shape ='s')
            plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"start.png")
            plt.close()
                     
    
            #########Setup Proposal
            ideal_population = sum(grid_partition["population"].values()) / len(grid_partition)
    
            tree_proposal = partial(recom,
                                   pop_col="population",
                                   pop_target=ideal_population,
                                   epsilon=0.05,
                                   node_repeats=1
                                  )
                                  
            #######BUILD MARKOV CHAINS
    
    
            exp_chain = MarkovChain(slow_reversible_propose_bi ,Validator([single_flip_contiguous,popbound#,boundary_condition
            ]), accept = always_accept, initial_state=grid_partition, 
            total_steps =  1000000)
    
    
    
    
    
            #########Run MARKOV CHAINS
    
            rsw = []
            rmm = []
            reg = []
            rce = []
            rbn=[]
            waits= []
    
            slopes = []
            angles = []
            lifts=[]
            loglifts=[]
    
            import time
    
            st = time.time()
    
    
            t=0
            for part in exp_chain:
                rce.append(len(part["cut_edges"]))
                #waits.append(part["geom"])
                rbn.append(len(list(part["b_nodes"])))
                #print(part["slope"])
                #temp = part["slope"]
                
                #enda = ((temp[0][0][0]+temp[0][1][0])/2,(temp[0][0][1]+temp[0][1][1])/2)  
                #endb = ((temp[1][0][0]+temp[1][1][0])/2,(temp[1][0][1]+temp[1][1][1])/2)  
                
                #if endb[0]!= enda[0]:
                
               #     slope = (endb[1]-enda[1])/(endb[0]-enda[0])
               # else:
               #     slope = np.Inf
                
               # slopes.append(slope)
                for edge in part["cut_edges"]:
                    graph[edge[0]][edge[1]]["cut_times"] += 1
                    #print(graph[edge[0]][edge[1]]["cut_times"])


    
                #anga = (temp[0][0]-20,temp[0][1]-20)
                #angb = (temp[1][0]-20,temp[1][1]-20)    
               # anga = (enda[0]-20,enda[1]-20)
                #angb = (endb[0]-20,endb[1]-20)    
               
               # angles.append(np.arccos(np.clip(np.dot(anga / np.linalg.norm(anga),angb / np.linalg.norm(angb)),-1,1)))
               
                if t%50000 == 0:
                    plt.figure()
                    nx.draw(graph, pos = {x:x for x in graph.nodes()}, node_color = [colors[dict(part.assignment)[x]] for x in graph.nodes()] ,node_size = ns, node_shape ='s')
                    plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"middle"+str(t)+".png")
                    plt.close()
                    
                    #for partname, subgraph in part.subgraphs.items():
                    L1 = nx.laplacian_matrix(part.subgraphs[1]).todense()
                    (sign1, logdet1) = np.linalg.slogdet(L1[1:,1:])
                    ctrees1 = np.exp(logdet1)#np.linalg.det(L1[1:,1:])
                    print(t,1,ctrees1)
                    print(t,1,logdet1)
                    Lm1 = nx.laplacian_matrix(part.subgraphs[-1]).todense()
                    (signm1, logdetm1) = np.linalg.slogdet(Lm1[1:,1:])
                    ctreesm1 = np.exp(logdetm1)#np.linalg.det(Lm1[1:,1:])
                    print(t,-1,ctreesm1)  
                    print(t,-1,logdetm1)        
                    lifts.append(ctrees1*ctreesm1*len(part["cut_edges"]))
                    loglifts.append(logdet1+logdetm1+np.log(len(part["cut_edges"])))
                    
                    print(t,'lifts',lifts[-1])
                    print(t,'loglifts',loglifts[-1])
                              
                    with open("./plots/States/gridtree/TREES"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"middle"+str(t)+".txt",'w') as treef:
                        treef.write(str(ctrees1))
                        treef.write("\n\n")
                        treef.write(str(ctreesm1))
                        treef.write("\n\n")
                        treef.write(str(len(part["cut_edges"])))
                        treef.write("\n\n")                        
                        treef.write(str(lifts[-1]))
                        treef.write("\n\n")
                        treef.write(str(logdet1))
                        treef.write("\n\n")
                        treef.write(str(logdetm1))
                        treef.write("\n\n")
                        treef.write(str(np.log(len(part["cut_edges"]))))
                        treef.write("\n\n")
                        treef.write("\n\n")
                        treef.write(str(loglifts[-1]))
                        
                        
                        
                        
                    np.savetxt("./plots/States/gridtree/LAP1"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"middle"+str(t)+".csv",L1, fmt='%d', delimiter = ',')
                    np.savetxt("./plots/States/gridtree/LAPm1"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"middle"+str(t)+".csv",Lm1, fmt='%d', delimiter = ',')
                    
                    
                    
                    
                
                if part.flips is not None:
                    f = list(part.flips.keys())[0]
                    graph.node[f]["part_sum"]=graph.node[f]["part_sum"]-part.assignment[f]*(t-graph.node[f]["last_flipped"])
                    graph.node[f]["last_flipped"]=t
                    graph.node[f]["num_flips"]=graph.node[f]["num_flips"]+1
    
                t+=1
                """
                plt.figure()
                nx.draw(graph, pos = {x:x for x in graph.nodes()}, node_color = [dict(part.assignment)[x] for x in graph.nodes()] ,node_size = ns, node_shape ='s',cmap = 'tab20')
                plt.savefig(f"./Figures/recom_{part['step_num']:02d}.png")
                plt.close()
                """
            print("finished no", st-time.time())
            with open("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"wait.txt",'w') as wfile:
                wfile.write(str(sum(waits)))

    


            for n in graph.nodes():
                if graph.node[n]["last_flipped"] == 0:
                    graph.node[n]["part_sum"]=t*part.assignment[n]
                graph.node[n]["lognum_flips"] = math.log(graph.node[n]["num_flips"] + 1) 

    
                        
                        
            #print(len(rsw[-1]))        
            #print(graph[(1,0)][(0,1)]["cut_times"])
            
            
            plt.figure()
            plt.plot(lifts,'ko')
            plt.ylabel("Tree Lifts")
            plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"lifts.png")
            plt.close()
            
            
            plt.figure()
            plt.plot(loglifts,'ko')
            plt.ylabel("Log Tree Lifts")
            plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"loglifts.png")
            plt.close()
            

            plt.figure()
            nx.draw(graph, pos = {x:x for x in graph.nodes()}, node_color = [0 for x in graph.nodes()] ,node_size = 1, edge_color = [graph[edge[0]][edge[1]]["cut_times"] for edge in graph.edges()], node_shape ='s',cmap = 'jet',width =5)
            plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"edges.png")
            plt.close()


    
            plt.figure()
            nx.draw(graph, pos = {x:x for x in graph.nodes()}, node_color = [colors[dict(part.assignment)[x]] for x in graph.nodes()] ,node_size = ns, node_shape ='s')
            plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"end.png")
            plt.close()
    
    
            A2 = np.zeros([50,50])
    
            for n in graph.nodes():
                A2[n[0],n[1]] = dict(part.assignment)[n]
                
                
            plt.figure()
            plt.imshow(A2,cmap='jet')
            plt.colorbar()
            plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"end2.png")
            plt.close()
    
    
            plt.figure()
            nx.draw(graph, pos = {x:x for x in graph.nodes()}, node_color = [graph.nodes[x]["part_sum"] for x in graph.nodes()] ,node_size = ns, node_shape ='s',cmap = 'jet')
            plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"wca.png")
            plt.close()
    
    
            A2 = np.zeros([50,50])
    
            for n in graph.nodes():
                A2[n[0],n[1]] = graph.nodes[n]["part_sum"]
                
                
            plt.figure()
            plt.imshow(A2,cmap='jet')
            plt.colorbar()

            plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"wca2.png")
            plt.close()
    
    
    

    
    
            plt.figure()
            plt.title("Flips")
            nx.draw(graph,pos= {x:x for x in graph.nodes()},node_color=[graph.nodes[x]["num_flips"] for x in graph.nodes()],node_size=ns,node_shape='s',cmap="jet")
            plt.title("Flips")
            plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"flip.png")
            plt.close()
    
    
            A2 = np.zeros([50,50])
    
            for n in graph.nodes():
                A2[n[0],n[1]] = graph.nodes[n]["num_flips"]
                
                
            plt.figure()
            plt.imshow(A2,cmap='jet')
            plt.colorbar()
            plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"flip2.png")
            plt.close()


            plt.figure()
            plt.title("Flips")
            nx.draw(graph,pos= {x:x for x in graph.nodes()},node_color=[graph.nodes[x]["lognum_flips"] for x in graph.nodes()],node_size=ns,node_shape='s',cmap="jet")
            plt.title("Flips")
            plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"logflip.png")
            plt.close()
    
    
            A2 = np.zeros([50,50])
    
            for n in graph.nodes():
                A2[n[0],n[1]] = graph.nodes[n]["lognum_flips"]
                
                
            plt.figure()
            plt.imshow(A2,cmap='jet')
            plt.colorbar()
            plt.savefig("./plots/States/gridtree/"+str(alignment)+"B"+str(int(100*base))+"P"+str(int(100*pop1))+"logflip2.png")
            plt.close()
