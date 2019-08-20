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
from gerrychain.tree import recursive_tree_part, bipartition_tree_random
#########################################################################
for plan_name in ['HB5005','REMEDIAL','HB7001','HB7002','HB7002ANS','HB7003','NAACP','PRINCETON']:
    run_param = plan_name

    newdir = "./Outputs/Plan_Stats/"+run_param+"/"
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

    df = gpd.read_file("./State_Data/VA_Trimmed_allplans.shp")
    
    graph.add_data(df, [plan_name])

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


    initial_partition = Partition(graph, plan_name, updater)
    print("built initial partition")


    ideal_population = sum(initial_partition["population"].values()) / len(
        initial_partition
    )
    # print(ideal_population)

    tree_proposal = partial(
        recom, pop_col="TOT_POP", pop_target=ideal_population, epsilon=0.01, node_repeats=1, method = bipartition_tree_random
    )

    compactness_bound = constraints.UpperBound(
        lambda p: len(p["cut_edges"]), 2 * len(initial_partition["cut_edges"])
    )



    print("initialized chain")
    
    #print([x for x in initial_partition.parts.keys()])
    #print(initial_partition["BVAP"].percents("First"))


    with open(newdir + "Start_Values.txt", "w") as f:
        f.write(f"Values for Starting Plan: {plan_name}\n \n ")
        f.write("Initial Cut: " + str(len(initial_partition["cut_edges"])))
        f.write("\n")
        f.write("\n")
        
        temp1 = [x for x in initial_partition.parts.keys()]
        temp2 = initial_partition["BVAP"].percents("First")
        temp3=[]
        for ind in range(len(temp1)):
            if temp1[ind] in  ['80', '89', '100', '79', '90', '83', '95', '94', '93', '97', '96', '64', '75', '62', '63', '91', '92', '70', '73', '74', '68', '72', '71', '66', '78', '76', '77', '81', '55', '61', '85', '69', '27']:
                temp3.append(temp2[ind])
                    

        for elect in range(num_elections):
        
        
            
                    
            f.write("RIGHT SET OF BVAPS: "+str(sorted(temp3)))
            
            f.write("\n")
            f.write("\n")
            

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

