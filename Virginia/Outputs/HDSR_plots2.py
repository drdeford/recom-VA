# sns.set_style('white')
import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


# sns.set_style('darkgrid')
sns.set_style("darkgrid", {"axes.facecolor": ".97"})


def draw_plot(data, offset, edge_color, fill_color):
    pos = 10*np.arange(data.shape[1])+offset
    #bp = ax.boxplot(data, positions= pos, widths=0.3, patch_artist=True, manage_xticks=False)
    bp = ax.boxplot(data, positions= pos,widths=.5, whis=[1,99],showfliers=False, patch_artist=True, manage_ticks=False,zorder=4)
    for element in ['boxes', 'whiskers', 'medians', 'caps']:
        plt.setp(bp[element], color=edge_color,zorder=4)
    for patch in bp['boxes']:
        patch.set(facecolor=fill_color,zorder=0)



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

newdir = "./Plots/Compare/"

datadir1= "./ReCOM_Enacted_uu/"
datadir2= "./ReCOM_Tree31_uutk3/"
datadir3= "./ReCOM_Tree99_uutk3/"

datadir1= "./FLIP_Enacted/"
datadir2= "./FLIP_Tree31/"
datadir3= "./FLIP_Tree99/"


os.makedirs(os.path.dirname(newdir + "init.txt"), exist_ok=True)
with open(newdir + "init.txt", "w") as f:
    f.write("Created Folder")

max_steps = 10000000#20000#
step_size = 10000#2000#

ts = [x * step_size for x in range(1, int(max_steps / step_size) + 1)]

for elect in range(1):
    a = []
    b = []
    c = []
    for t in ts:
        tempvotes = np.loadtxt(
            datadir1 + election_names[elect] + "_" + str(t) + ".csv", delimiter=","
        )
        for s in range(step_size):
            a.append(tempvotes[s, :])

        tempvotes = np.loadtxt(
            datadir2 + election_names[elect] + "_" + str(t) + ".csv", delimiter=","
        )
        for s in range(step_size):
            b.append(tempvotes[s, :])

        tempvotes = np.loadtxt(
            datadir3 + election_names[elect] + "_" + str(t) + ".csv", delimiter=","
        )
        for s in range(step_size):
            c.append(tempvotes[s, :])

    a = np.array(a)
    b = np.array(b)
    c = np.array(c)

    #medianprops = dict(color="black")

    fig, ax = plt.subplots()
    draw_plot(a,-2,'r',None)
    draw_plot(b,0,'y',None)
    draw_plot(c,2,'b',None)

    plt.plot([],[],color='r',label='Enacted')
    plt.plot([],[],color='y',label='Seed 31')
    plt.plot([],[],color='b',label='Seed 99')

    plt.legend()

    #plt.xticks([5,10,15,20,25,30],[5,10,15,20,25,30])
    plt.xticks([],[])
    #plt.xlim([.5,34])
    plt.xlabel("Indexed Districts")
    plt.ylabel("BVAP %")

    plt.legend()

    plt.savefig("./Plots/HDSR2/FLIPseed_compare.png")
    fig = plt.gcf()
    fig.set_size_inches((12,8), forward=False)
    fig.savefig("./Plots/HDSR2/FLIPseed_compare2.png", dpi=1000)
    plt.close()



