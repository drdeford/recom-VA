import geopandas as gpd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import json
import os

hb5005=[0.07253105784916551, 0.10824003731840597, 0.13398271190814087, 0.13546485835604286, 0.13711169698855355, 0.15123053901747907, 0.16064167834079338, 0.1671916890080429, 0.1714120178778348, 0.18436445635035012, 0.18601774940250362, 0.18926480993117642, 0.19611679964912873, 0.2102159841056207, 0.2258110938774205, 0.24239331194711772, 0.2456163350922618, 0.2514435871257134, 0.27596109603820584, 0.2945835079944907, 0.33525600505689, 0.5518884518212926, 0.5534953948361769, 0.5542911182914335, 0.5545930898968396, 0.5629610159189105, 0.5636955706345689, 0.5658926317188226, 0.5723899599854493, 0.5877664632354213, 0.5952869519900984, 0.5997255028212211, 0.6071539251985842]

#6110

hb7001=[0.10824003731840597, 0.13539411689800077, 0.13711169698855355, 0.15831392290930948, 0.15964176113233544, 0.1659980607258969, 0.1714120178778348, 0.17896382948475803, 0.18391390540923253, 0.20830404835535563, 0.2112792472024415, 0.2258110938774205, 0.22767239760290928, 0.24239331194711772, 0.25487725928647664, 0.2732269699021939, 0.2802270884022709, 0.3047803894000061, 0.3342634208283905, 0.33525600505689, 0.4240456445754232, 0.4797818705628151, 0.4936303619414126, 0.5012558819788885, 0.5085824894818511, 0.514421058943455, 0.5180834320390713, 0.5225580889610064, 0.5269443588037711, 0.5294892751109184, 0.5353813664186369, 0.5542911182914335, 0.5726470920201874]

#6023

hb7002 = [0.11272215546138764, 0.11382279273638071, 0.13711169698855355, 0.1431238231283683, 0.14966209566301425, 0.15226036331244233, 0.1671916890080429, 0.16792228057517744, 0.18164005929156798, 0.18498686947921286, 0.19566728639560263, 0.19776340017230604, 0.204824382967182, 0.2223996902376458, 0.2546487352527067, 0.27262543679256324, 0.27400258839695146, 0.27636154803305363, 0.33525600505689, 0.34072008366117223, 0.38048553478631814, 0.4427302631578947, 0.46755868865922934, 0.47420802913877663, 0.5144438010680049, 0.5337139113752282, 0.5440908795423519, 0.5509299412161511, 0.5526805233791235, 0.5644361739448758, 0.5859314139857126, 0.6002399960000666, 0.6176769387064003]

#5300

hb7002ans=[0.11272215546138764, 0.11382279273638071, 0.13711169698855355, 0.1431238231283683, 0.14966209566301425, 0.15226036331244233, 0.1671916890080429, 0.16792228057517744, 0.1821179670236274, 0.18498686947921286, 0.19566728639560263, 0.19776340017230604, 0.204824382967182, 0.23055644901897715, 0.2546487352527067, 0.27262543679256324, 0.27636154803305363, 0.2918468041336526, 0.33525600505689, 0.34072008366117223, 0.38048553478631814, 0.4427302631578947, 0.46755868865922934, 0.47420802913877663, 0.5144438010680049, 0.5337139113752282, 0.5440908795423519, 0.5509299412161511, 0.5526805233791235, 0.5644361739448758, 0.5730348308430363, 0.5859314139857126, 0.6176769387064003]

#5294

hb7003 = [0.13572773958819648, 0.13711169698855355, 0.14796551345069647, 0.15095212906638455, 0.15291473804253036, 0.15540199191685913, 0.15880834295808138, 0.16379254351779682, 0.17304222905251904, 0.1933225172914354, 0.2056550669464916, 0.2097701631154549, 0.2102159841056207, 0.23731004953404417, 0.2618367011664202, 0.27242391884676986, 0.2777513667696696, 0.28192966219342896, 0.29736540294312436, 0.33525600505689, 0.46796638988394934, 0.4741456492727672, 0.483924460257917, 0.4936303619414126, 0.514421058943455, 0.5180834320390713, 0.5224328799142941, 0.5269443588037711, 0.5294892751109184, 0.5353813664186369, 0.5451518721507965, 0.5542911182914335, 0.5726470920201874]

#6008

naacp = [0.07770614801864802, 0.10824003731840597, 0.13711169698855355, 0.1589873417721519, 0.16230064689681675, 0.16532213850467584, 0.1671916890080429, 0.16999477943095798, 0.17421161860731943, 0.18237143385505103, 0.22424945810790053, 0.2258110938774205, 0.2483963688919924, 0.27596109603820584, 0.27736367946894264, 0.2843588759541171, 0.2933731186674099, 0.32515791700667934, 0.33525600505689, 0.42489150904773604, 0.4252847965738758, 0.43831846964824456, 0.4596767617463346, 0.4649128675117905, 0.4750695811553866, 0.49721274921301156, 0.49951068840025686, 0.5036721385566478, 0.506660047985439, 0.5212373624976155, 0.5459902316054829, 0.5759931611800191, 0.5852229871776354]

#4874

princeton = [0.11285994764397905, 0.11923800918605526, 0.13897052649256864, 0.15490575550098393, 0.16074647274535397, 0.1644754539340955, 0.1924158221641059, 0.20748078731476377, 0.2389974186557889, 0.247833944163894, 0.2492809364548495, 0.24964654261802227, 0.26311193301484564, 0.26698328522521314, 0.3289734116249918, 0.32909409158231023, 0.3425087866629343, 0.3499765322802531, 0.3698742706878133, 0.38867736303932165, 0.40328054298642535, 0.4191690178691388, 0.43597163215468504, 0.4438127358869349, 0.4490251331585034, 0.451499388004896, 0.463070511068598, 0.46434859154929575, 0.4681934272618101, 0.49059480627868907, 0.5000960491771788, 0.5379312500977005, 0.570986567004151]

#4559

remedial=[0.07253105784916551, 0.10824003731840597, 0.13546485835604286, 0.13711169698855355, 0.1538334625728003, 0.1671916890080429, 0.1685010703028873, 0.18184318229613214, 0.18436445635035012, 0.2129347948556988, 0.23096673364590312, 0.25339647673266497, 0.272224351390532, 0.27379063261733766, 0.27596109603820584, 0.31127146056232896, 0.3145810997089907, 0.3231048121292024, 0.32515968477810037, 0.33525600505689, 0.4022752330359828, 0.4193422274348825, 0.42884963613379345, 0.47364768859663614, 0.47474830423571535, 0.5137779093091316, 0.522929968159768, 0.5244954056388911, 0.5386911037223259, 0.5401329975722665, 0.5436851282902931, 0.5437768858042953, 0.5492214644189937]



#4708


tree31=[0.09372128134436129, 0.135168609835849, 0.14111654377662836, 0.16053649304698428, 0.1680385230605172, 0.16936988904346645, 0.17919346982650028, 0.18052098078456258, 0.19916656864886664, 0.21587935652685525, 0.223407917383821, 0.26329642982484125, 0.2688731258537549, 0.2763703213945913, 0.2777111698720581, 0.2865828425911189, 0.30190981342461576, 0.3094254357650097, 0.3216829116208854, 0.32765688006131904, 0.32804753876713016, 0.3592459834858862, 0.37433767959471803, 0.4034688314411527, 0.44000389111366917, 0.45793038180581436, 0.4643128970627068, 0.4973517039358278, 0.5237142715445273, 0.5770066146707473, 0.631588696911522, 0.6827740566904408, 0.6876783313108786]

#6558

tree99=[0.0897465026453309, 0.12868573934012933, 0.12936187264166696, 0.13786807812133597, 0.14722680120552142, 0.17933390264731, 0.20061398599026442, 0.2231551465729803, 0.22575038223223626, 0.22716061482171152, 0.24196774193548387, 0.2587116900013306, 0.2790637317718049, 0.28488719980593513, 0.285673952036139, 0.2900348055283344, 0.30937338512060547, 0.3733133942833837, 0.37990591176760863, 0.3860371546223242, 0.3881455418915033, 0.389077225024736, 0.40379168276039656, 0.4104842801263482, 0.41143578028572403, 0.41255241444609764, 0.48576306019623444, 0.49806962712884195, 0.5044071780499932, 0.5405626020344091, 0.5481079519000754, 0.5678214818579791, 0.5716793354752763]

#6603

vec =[hb5005,hb7001,hb7002,hb7002ans,hb7003,naacp,princeton,remedial]

for l in vec:
    print(len(l))
    
def draw_plot(data, offset, edge_color, fill_color):
    pos = np.arange(data.shape[1])+offset
    #bp = ax.boxplot(data, positions= pos, widths=0.3, patch_artist=True, manage_xticks=False)
    bp = ax.boxplot(data, positions= pos,widths=.5, whis=[1,99],showfliers=False, patch_artist=True, manage_ticks=False,zorder=4)
    for element in ['boxes', 'whiskers', 'medians', 'caps']:
        plt.setp(bp[element], color=edge_color,zorder=4)
    for patch in bp['boxes']:
        patch.set(facecolor=fill_color,zorder=0)

datadir="./ReCOM_Enacted_uu/"
max_steps = 20000#10000000#
step_size = 2000#10000#

ts = [x * step_size for x in range(1, int(max_steps / step_size) + 1)]
size=0
for elect in range(1):
    a = []
    for t in ts:
        tempvotes = np.loadtxt(
            datadir + "BVAP_" + str(t) + ".csv", delimiter=","
        )
        for s in range(step_size):
            if tempvotes[s,-1]<.6:
                a.append(tempvotes[s,:])
                size+=1

    a = np.array(a)
    print(size)
    
    
fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.add_patch(patches.Rectangle((0, .37), 35, .18,color='honeydew'))
plt.plot([0,34], [.55, .55], 'lightgreen')
plt.plot([0,34], [.37, .37], 'lightgreen')
plt.axvline(x=4.5, color='gray', linestyle='--')
plt.axvline(x=8.5, color='gray', linestyle='--')
plt.axvline(x=17.5, color='gray', linestyle='--')
plt.axvline(x=21.5, color='gray', linestyle='--')



medianprops = dict(color="black")

plt.boxplot(a,whis=[1,99],showfliers=False,medianprops=medianprops)


plt.plot(range(1,34),hb5005,'o',color='r', label = 'Enacted')
plt.plot(range(1,34),hb7002,'o',color='g', label = 'GOP1')
plt.plot(range(1,34),hb7002ans,'o',color='y', label = 'GOP2')
plt.plot(range(1,34),hb7003,'o',color='purple', label = 'GOP3')

#plt.xticks(range(1,34),range(1,34))
plt.xticks([5,10,15,20,25,30],[5,10,15,20,25,30])
plt.xlim([.5,34])
plt.xlabel("Indexed Districts")
plt.ylabel("BVAP %")
plt.legend()

plt.savefig("./Plots/HDSR2/REBP601.png")
fig = plt.gcf()
fig.set_size_inches((12,8), forward=False)
fig.savefig("./Plots/HDSR2/REBP6012.png", dpi=1000)
plt.close()

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.add_patch(patches.Rectangle((0, .37), 35, .18,color='honeydew'))
plt.plot([0,34], [.55, .55], 'lightgreen')
plt.plot([0,34], [.37, .37], 'lightgreen')
plt.axvline(x=4.5, color='gray', linestyle='--')
plt.axvline(x=8.5, color='gray', linestyle='--')
plt.axvline(x=17.5, color='gray', linestyle='--')
plt.axvline(x=21.5, color='gray', linestyle='--')



medianprops = dict(color="black")

plt.boxplot(a,whis=[1,99],showfliers=False,medianprops=medianprops)


plt.plot(range(1,34),hb7001,'o',color='b', label = 'DEM')
plt.plot(range(1,34),naacp,'o',color='orange', label = 'NAACP')
plt.plot(range(1,34),princeton,'o',color='tan', label = 'Princeton')
plt.plot(range(1,34),remedial,'o',color='hotpink', label = 'Remedial')

#plt.xticks(range(1,34),range(1,34))
plt.xticks([5,10,15,20,25,30],[5,10,15,20,25,30])
plt.xlim([.5,34])
plt.xlabel("Indexed Districts")
plt.ylabel("BVAP %")
plt.legend()

plt.savefig("./Plots/HDSR2/REBP602.png")
fig = plt.gcf()
fig.set_size_inches((12,8), forward=False)
fig.savefig("./Plots/HDSR2/REBP6022.png", dpi=1000)
plt.close()

