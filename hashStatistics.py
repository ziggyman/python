import csvFree,csvData
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
#cs = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/hash_CSCoords_051021.csv',',',False)

print('hello')
nESO2000 = 1510
nESOT = 1268
nESOL = 42
nESOP = 22
nESO = 1332
nK = 1510
nHASH = 4631
nHASHgPNe = 3804
nHASHT = 3299
nHASHL = 455
nHASHP = 877
nDSH = 119+16
nDSHT = 86
nDSHL = 10+6
nDSHP = 23+10
nFrenchAmateursT = 110
nFrenchAmateursL = 36
nFrenchAmateursP = 32
nFrenchAmateurs = nFrenchAmateursT+nFrenchAmateursL+nFrenchAmateursP
nIPHAST = 119
nIPHASL = 18
nIPHASP = 16
nIPHAS = nIPHAST+nIPHASL+nIPHASP
nOthers = 970
nMASH = 1190

#plt.plot([0,4],[0,0])
eso = PatchCollection([Rectangle((0.,0.),nESO,1.)],
                      facecolor='g',
                      alpha=1.,
                      edgecolor=None)
koh = PatchCollection([Rectangle((nESO,0.),nK-nESO,1.)],
                      facecolor='b',
                      alpha=1.,
                      edgecolor=None)
ama = PatchCollection([Rectangle((nK,0.),nFrenchAmateurs+nDSH,1.)],
                      facecolor='m',
                      alpha=1.,
                      edgecolor=None)
iph = PatchCollection([Rectangle((nK+nFrenchAmateurs+nDSH,0),nIPHAS,1.)],
                      facecolor='r',
                      alpha=1.,
                      edgecolor=None)
oth = PatchCollection([Rectangle((nK+nFrenchAmateurs+nDSH+nIPHAS,0),nOthers-nFrenchAmateurs-nDSH-nIPHAS,1.)],
                        facecolor='y',
                        alpha=0.9,
                        edgecolor=None)
mas = PatchCollection([Rectangle((nK+nOthers,0),nMASH,1.)],
                        facecolor='orange',
                        alpha=0.9,
                        edgecolor=None)
has = PatchCollection([Rectangle((nK+nOthers+nMASH,0),nHASH-nK-nOthers-nMASH,1.)],
                        facecolor='c',
                        alpha=0.9,
                        edgecolor=None)
hasg = PatchCollection([Rectangle((nK+nOthers+nMASH,0),nHASHgPNe-nK-nOthers-nMASH,1.)],
                        facecolor='c',
                        alpha=0.9,
                        edgecolor=None)
fig, ax = plt.subplots(1,figsize=(nHASH/3.,10.),dpi=40)
ax.add_collection(eso)
ax.add_collection(koh)
ax.add_collection(ama)
ax.add_collection(iph)
ax.add_collection(oth)
ax.add_collection(mas)
ax.add_collection(has)
ax.add_collection(has)
plt.xlim([0, nHASH])
plt.ylim(0,1.)
plt.yticks([])
#plt.yaxis('off')
plt.xticks(fontsize=64)
plt.savefig('/Users/azuri/daten/uni/HKU/APN8e/nPNe.png', bbox_inches='tight', transparent=True)
plt.show()
plt.close()
eso_a = PatchCollection([Rectangle((0.,0.),nESO,1.)],
                      facecolor='g',
                      alpha=1.,
                      edgecolor=None)
koh_a = PatchCollection([Rectangle((nESO,0.),nK-nESO,1.)],
                      facecolor='b',
                      alpha=1.,
                      edgecolor=None)
ama_a = PatchCollection([Rectangle((nK,0.),nFrenchAmateurs+nDSH,1.)],
                      facecolor='m',
                      alpha=1.,
                      edgecolor=None)
iph_a = PatchCollection([Rectangle((nK+nFrenchAmateurs+nDSH,0),nIPHAS,1.)],
                      facecolor='r',
                      alpha=1.,
                      edgecolor=None)
oth_a = PatchCollection([Rectangle((nK+nFrenchAmateurs+nDSH+nIPHAS,0),nOthers-nFrenchAmateurs-nDSH-nIPHAS,1.)],
                        facecolor='y',
                        alpha=0.9,
                        edgecolor=None)
mas_a = PatchCollection([Rectangle((nK+nOthers,0),nMASH,1.)],
                        facecolor='orange',
                        alpha=0.9,
                        edgecolor=None)
hasg_a = PatchCollection([Rectangle((nK+nOthers+nMASH,0),nHASHgPNe-nK-nOthers-nMASH,1.)],
                        facecolor='c',
                        alpha=0.9,
                        edgecolor=None)

fig, ax = plt.subplots(1,figsize=(nHASHgPNe/3.,8.),dpi=40)
ax.add_collection(eso_a)
ax.add_collection(koh_a)
ax.add_collection(ama_a)
ax.add_collection(iph_a)
ax.add_collection(oth_a)
ax.add_collection(mas_a)
ax.add_collection(hasg_a)
plt.xlim([0, nHASH])
plt.ylim(0,1.)
plt.yticks([])
#plt.yaxis('off')
plt.xticks(fontsize=64)
plt.savefig('/Users/azuri/daten/uni/HKU/APN8e/ngPNe.png', bbox_inches='tight', transparent=True)
plt.show()

labels = 'ESO','Kohoutek','Amateurs','IPHAS','MASH','New HASH','Others'
explode = (0, 0, 0, 0.01, 0.01, 0.01, 0)
fig1,ax1 = plt.subplots()
ax1.pie([nESO,nK-nESO,nFrenchAmateurs+nDSH,nIPHAS,nMASH,nHASHgPNe-nK-nOthers-nMASH,nOthers-nFrenchAmateurs-nDSH-nIPHAS],
        labels=labels,
        autopct='%1.1f%%',
        shadow=False,
        #explode=explode, 
        startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.savefig('/Users/azuri/daten/uni/HKU/APN8e/catalogues_gPNe.png', bbox_inches='tight', transparent=True)
plt.show()

hashCSVFile = '/Users/azuri/daten/uni/HKU/APN8e/hash_tlp.csv'

hash = csvFree.readCSVFile(hashCSVFile)

nR = len(hash.find('mainClass','R'))
nE = len(hash.find('mainClass','E'))
nB = len(hash.find('mainClass','B'))
nI = len(hash.find('mainClass','I'))
nA = len(hash.find('mainClass','A'))
nS = len(hash.find('mainClass','S'))

labels = 'R','E','B','I','A','S'

fig1,ax1 = plt.subplots()
ax1.pie([nR,nE,nB,nI,nA,nS], labels=labels, autopct='%1.1f%%',
        shadow=False, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.savefig('/Users/azuri/daten/uni/HKU/APN8e/morphologies_gPNe.png', bbox_inches='tight', transparent=True)
plt.show()


nT = len(hash.find('PNstat','T'))
nL = len(hash.find('PNstat','L'))
nP = len(hash.find('PNstat','P'))

labels = 'T','L','P'

fig1,ax1 = plt.subplots()
ax1.pie([nT,nL,nP], labels=labels, autopct='%1.1f%%',
        shadow=False, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.savefig('/Users/azuri/daten/uni/HKU/APN8e/PNstat_gPNE.png', bbox_inches='tight', transparent=True)
plt.show()
