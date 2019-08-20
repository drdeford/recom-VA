import csv
import os
from functools import partial
import json
import numpy as np
import geopandas as gpd
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from networkx.readwrite import json_graph
import random

from gerrychain import (
    Election,
    Graph,
    MarkovChain,
    Partition,
    accept,
    constraints,
    updaters,
)
from gerrychain.constraints import single_flip_contiguous
from gerrychain.metrics import efficiency_gap, mean_median
from gerrychain.proposals import recom, propose_random_flip
from gerrychain.updaters import cut_edges
from gerrychain.tree import recursive_tree_part, bipartition_tree_random, PopulatedGraph,contract_leaves_until_balanced_or_none,find_balanced_edge_cuts
#########################################################################
def get_spanning_tree_u_w(G):
    node_set=set(G.nodes())
    x0=random.choice(tuple(node_set))
    x1=x0
    while x1==x0:
        x1=random.choice(tuple(node_set))
    node_set.remove(x1)
    tnodes ={x1}
    tedges=[]
    current=x0
    current_path=[x0]
    current_edges=[]
    while node_set != set():
        next=random.choice(list(G.neighbors(current)))
        current_edges.append((current,next))
        current = next
        current_path.append(next)

        if next in tnodes:
            for x in current_path[:-1]:
                node_set.remove(x)
                tnodes.add(x)
            for ed in current_edges:
                tedges.append(ed)
            current_edges = []
            if node_set != set():
                current=random.choice(tuple(node_set))
            current_path=[current]


        if next in current_path[:-1]:
            current_path.pop()
            current_edges.pop()
            for i in range(len(current_path)):
                if current_edges !=[]:
                    current_edges.pop()
                if current_path.pop() == next:
                    break
            if len(current_path)>0:
                current=current_path[-1]
            else:
                current=random.choice(tuple(node_set))
                current_path=[current]

    #tgraph = Graph()
    #tgraph.add_edges_from(tedges)
    return G.edge_subgraph(tedges)

def my_uu_bipartition_tree_random(
    graph,
    pop_col,
    pop_target,
    epsilon,
    node_repeats=1,
    spanning_tree=None,
    choice=random.choice):
    populations = {node: graph.nodes[node][pop_col] for node in graph}

    possible_cuts = []
    if spanning_tree is None:
        spanning_tree = get_spanning_tree_u_w(graph)

    while len(possible_cuts) == 0:
        spanning_tree = get_spanning_tree_u_w(graph)
        h = PopulatedGraph(spanning_tree, populations, pop_target, epsilon)
        possible_cuts = find_balanced_edge_cuts(h, choice=choice)

    return choice(possible_cuts).subset


run_param = 'ReCOM_Enacted_uu'

newdir = "./Outputs/"+run_param+"/"
os.makedirs(os.path.dirname(newdir + "init.txt"), exist_ok=True)
with open(newdir + "init.txt", "w") as f:
    f.write("Created Folder")


def b_nodes(partition):
    return {(x[0], partition.assignment[x[1]]) for x in partition["cut_edges"]
               }.union({(x[1], partition.assignment[x[0]]) for x in partition["cut_edges"]})

def slow_reversible_propose(partition):

    flip = random.choice(list(partition["b_nodes"]))
    
    return partition.flip({flip[0]: flip[1]})

#graph = Graph.from_file("./State_Data/VA_Trimmed.shp")


#with open("./State_Data/VA_Trimmedj.json", 'w') as outfile:
#    json.dump(json_graph.adjacency_data(graph), outfile)

#graph.to_json("./State_Data/VA_Trimmed.json")

graph = Graph.from_json("./State_Data/VA_Trimmed.json")

df = gpd.read_file("./State_Data/VA_Trimmed.shp")

print("loaded data")
#print("nodes",len(graph.nodes()))
#print("edges",len(graph.edges()))


for node in graph.nodes():
    graph.nodes[node]["nBVAP"] = graph.nodes[node]["VAPERSONS"] - graph.nodes[node]["VAPBLACK"]

"""
updater = { "population": updaters.Tally("TAPERSONS", alias="population"),
"cut_edges",cut_edges,
"BVAP":Election("BVAP",{"BVAP":"VAPBLACK","nBVAP":"nBVAP"}),
"LTGOV":Election("LTGOV",{"D":"D_LTGOV","R":"R_LTGOV"}),
"GOV":Election("GOV",{"D":"D_GOV","R":"R_GOV"}),
"AG":Election("AG",{"D":"D_ATTGEN","R":"R_ATTGEN"}),
}
"""



#print("parts",len(initial_partition))
#print(sorted(initial_partition["BVAP"].percents("BVAP")))

compactness_bound = constraints.UpperBound(
    lambda p: len(p["cut_edges"]), 12000)




num_elections = 4


election_names = [
    "BVAP",
    "LTGOV",
    "GOV",
    "AG"]

election_columns = [
    ["VAPBLACK", "nBVAP"],
    ["D_LTGOV", "R_LTGOV"],
    ["D_GOV", "R_GOV"],
    ["D_ATTGEN", "R_ATTGEN"],
    ["PRES12D", "PRES12R"],
    ]


updater = {
    "population": updaters.Tally("TAPERSONS", alias="population"),
    "cut_edges": cut_edges,
    'b_nodes':b_nodes,
}

elections = [
    Election(
        election_names[i],
        {"First": election_columns[i][0], "Second": election_columns[i][1]},
    )
    for i in range(num_elections)
]

election_updaters = {election.name: election for election in elections}

updater.update(election_updaters)


initial_partition = Partition(graph, "HB5005", updater)
print("built initial partition")


ideal_population = sum(initial_partition["population"].values()) / len(
    initial_partition
)
# print(ideal_population)

tree_proposal = partial(
    recom, pop_col="TAPERSONS", pop_target=ideal_population, epsilon=0.005, node_repeats=1, method =my_uu_bipartition_tree_random)# bipartition_tree_random
#)

compactness_bound = constraints.UpperBound(
    lambda p: len(p["cut_edges"]), 2 * len(initial_partition["cut_edges"])
)

chain = MarkovChain(
    proposal=tree_proposal,
    constraints=[
        #constraints.within_percent_of_ideal_population(initial_partition, 0.013),
        compactness_bound#,  single_flip_contiguous#no_more_discontiguous
    ],
    accept=accept.always_accept,
    initial_state=initial_partition,
    total_steps=20000,
)

print("initialized chain")


with open(newdir + "Start_Values.txt", "w") as f:
    f.write("Values for Starting Plan: Tree Recursive\n \n ")
    f.write("Initial Cut: " + str(len(initial_partition["cut_edges"])))
    f.write("\n")
    f.write("\n")

    for elect in range(num_elections):
        f.write(
            election_names[elect]
            + "District Percentages"
            + str(
                sorted(initial_partition[election_names[elect]].percents("First"))
            )
        )
        f.write("\n")
        f.write("\n")

        f.write(
            election_names[elect]
            + "Mean-Median :"
            + str(mean_median(initial_partition[election_names[elect]]))
        )

        f.write("\n")
        f.write("\n")

        f.write(
            election_names[elect]
            + "Efficiency Gap :"
            + str(efficiency_gap(initial_partition[election_names[elect]]))
        )

        f.write("\n")
        f.write("\n")

        f.write(
            election_names[elect]
            + "How Many Seats :"
            + str(initial_partition[election_names[elect]].wins("First"))
        )

        f.write("\n")
        f.write("\n")


print("wrote starting values")

pop_vec = []
cut_vec = []
votes = [[], [], [], []]
mms = []
egs = []
hmss = []


t = 0
for part in chain:

    pop_vec.append(sorted(list(part["population"].values())))
    cut_vec.append(len(part["cut_edges"]))
    mms.append([])
    egs.append([])
    hmss.append([])

    for elect in range(num_elections):
        votes[elect].append(sorted(part[election_names[elect]].percents("First")))
        mms[-1].append(mean_median(part[election_names[elect]]))
        egs[-1].append(efficiency_gap(part[election_names[elect]]))
        hmss[-1].append(part[election_names[elect]].wins("First"))

    t += 1
    if t%100 ==0:
        print(t)
    if t % 2000 == 0:
        print(t)
        with open(newdir + "mms" + str(t) + ".csv", "w") as tf1:
            writer = csv.writer(tf1, lineterminator="\n")
            writer.writerows(mms)

        with open(newdir + "egs" + str(t) + ".csv", "w") as tf1:
            writer = csv.writer(tf1, lineterminator="\n")
            writer.writerows(egs)

        with open(newdir + "hmss" + str(t) + ".csv", "w") as tf1:
            writer = csv.writer(tf1, lineterminator="\n")
            writer.writerows(hmss)

        with open(newdir + "pop" + str(t) + ".csv", "w") as tf1:
            writer = csv.writer(tf1, lineterminator="\n")
            writer.writerows(pop_vec)

        with open(newdir + "cuts" + str(t) + ".csv", "w") as tf1:
            writer = csv.writer(tf1, lineterminator="\n")
            writer.writerows([cut_vec])

        with open(newdir + "assignment" + str(t) + ".json", "w") as jf1:
            json.dump(dict(part.assignment), jf1)

        for elect in range(num_elections):
            with open(
                newdir + election_names[elect] + "_" + str(t) + ".csv", "w"
            ) as tf1:
                writer = csv.writer(tf1, lineterminator="\n")
                writer.writerows(votes[elect])

        df["plot" + str(t)] = df.index.map(dict(part.assignment))
        df.plot(column="plot" + str(t), cmap="tab20")
        plt.savefig(newdir + "plot" + str(t) + ".png")
        plt.close()

        votes = [[], [], [], []]
        mms = []
        egs = []
        hmss = []
        pop_vec = []
        cut_vec = []





