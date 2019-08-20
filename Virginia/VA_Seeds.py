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


from gerrychain import (
    Election,
    Graph,
    MarkovChain,
    Partition,
    accept,
    constraints,
    updaters,
)
from gerrychain.metrics import efficiency_gap, mean_median
from gerrychain.proposals import recom, propose_random_flip
from gerrychain.updaters import cut_edges
from gerrychain.tree import recursive_tree_part

run_param = 'build'

newdir = "./Outputs/"+run_param+"/"
os.makedirs(os.path.dirname(newdir + "init.txt"), exist_ok=True)
with open(newdir + "init.txt", "w") as f:
    f.write("Created Folder")


#graph = Graph.from_file("./State_Data/VA_Trimmed.shp")


#with open("./State_Data/VA_Trimmedj.json", 'w') as outfile:
#    json.dump(json_graph.adjacency_data(graph), outfile)

#graph.to_json("./State_Data/VA_Trimmed.json")

graph = Graph.from_json("./State_Data/VA_Trimmed.json")

df = gpd.read_file("./State_Data/VA_Trimmed.shp")

print("nodes",len(graph.nodes()))
print("edges",len(graph.edges()))

num_districts =33

for i in range(100):
    cddict =  recursive_tree_part(graph,range(num_districts),df["TAPERSONS"].sum()/num_districts,"TAPERSONS", .001,1)
    
    df["initial"]=df.index.map(cddict)
    
    df.plot(column="initial",cmap="jet")
    plt.savefig(newdir+"initial.png")
    plt.close()

    with open(newdir+"init"+str(i)+".json", 'w') as jf1:
        json.dump(cddict, jf1)

    print("saved Initial",i)

