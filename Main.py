import matplotlib.pyplot as plt
# import seaborn as sns
import numpy as np
import collections

# with open ("atomValuesShort.txt") as f:
with open ("atomValuesPerAA.txt") as f:
    contents = f.readlines()

database = {}
resultDB = {}
AA = ""
for line in contents:
    line = line.strip('\n')
    line = line.split("\t")
    results = []
    curAtom = ""
    if line[0].startswith("Amino Acid"):
        AA = line[0][len(line)-2]
        database[AA] = {}
        resultDB[AA] = {}

    else:
        results = []
        curAtom = line[0].strip()
        for j in range(1,len(line)-1):
            pair = line[j].split(",")
            results.append(pair)
        database[AA][curAtom] = results

for aminoAcid, atomTypes in database.items():

    for atom, values in atomTypes.items():
        sum = 0
        sumMinusCalpha = 0
        for pair in values:
            sum+= float(pair[1])
            ratio = (float(pair[1]) / float(pair[0])) if float(pair[0])!=0.0 else 0
            sumMinusCalpha += ratio
        avg = sum / len(values)
        avgMinusCalpha = sumMinusCalpha  / len(values)
        resultDB[aminoAcid][atom] = [avg,avgMinusCalpha]



acids = list(resultDB.keys())

fig, axarr = plt.subplots(10,2, sharey=True)
axlist = axarr.flatten()

i=0
for acid in acids:
    resultDB[acid] = collections.OrderedDict(sorted(resultDB[acid].items()))
    atoms = list(resultDB[acid].keys())
    values = list(resultDB[acid].values())
    intAvg = [itm[0] for itm in values]
    intNorm = [itm[1] for itm in values]
    N = len(intAvg)
    ind = np.arange(N)
    width = 0.35
    rects1 = axlist[i].bar(ind, intAvg, width, color='r')
    rects2 = axlist[i].bar(ind+width, intNorm, width, color='y')

    axlist[i].set_ylabel('Intensity Value', fontsize=40)
    axlist[i].set_xlabel('AA: '+ acid, fontsize=40)
    axlist[i].set_xticks(ind+width)
    axlist[i].set_xticklabels(atoms, fontsize=20)

    i+=1

plt.suptitle('Intensity Values Per Atom Type', fontsize=70)
fig.legend((rects1[0], rects2[0]), ('Avg Intensity', 'Atom / CA'), loc=1, fontsize=40, borderaxespad=8)
fig.set_figheight(70)
fig.set_figwidth(30)

plt.savefig("All.png")



### gln vs glu
fig, axarr = plt.subplots(1,2, sharey=True)

gluAtoms = list(resultDB['E'].keys())
glnAtoms = list(resultDB['Q'].keys())

gluValues = list(resultDB['E'].values())
glnValues = list(resultDB['Q'].values())

gluIntAvg = [itm[0] for itm in gluValues]
gluIntNorm = [itm[1] for itm in gluValues]

glnIntAvg = [itm[0] for itm in glnValues]
glnIntNorm = [itm[1] for itm in glnValues]

N = len(glnIntAvg)
ind = np.arange(N)
width = 0.35

rects1 = axarr[0].bar(ind, glnIntAvg, width, color='r')
rects2 = axarr[0].bar(ind+width, glnIntNorm, width, color='y')


axarr[0].set_ylabel('Intensity Value', fontsize=40)
axarr[0].set_xlabel('Glutamine', fontsize=40)
axarr[0].set_xticks(ind+width)
axarr[0].set_xticklabels(glnAtoms, fontsize=20)

N = len(gluIntAvg)
ind = np.arange(N)

rects1 = axarr[1].bar(ind, gluIntAvg, width, color='r')
rects2 = axarr[1].bar(ind+width, gluIntNorm, width, color='y')

axarr[1].set_ylabel('Intensity Value', fontsize=40)
axarr[1].set_xlabel('Glutamic Acid', fontsize=40)
axarr[1].set_xticks(ind+width)
axarr[1].set_xticklabels(gluAtoms, fontsize=20)




fig.set_figheight(20)
fig.set_figwidth(40)
plt.savefig("GluVsGln.png")



### asp vs asn
fig, axarr = plt.subplots(1,2, sharey=True)
asnAtoms = list(resultDB['N'].keys())
aspAtoms = list(resultDB['D'].keys())

asnValues = list(resultDB['N'].values())
aspValues = list(resultDB['D'].values())

asnIntAvg = [itm[0] for itm in asnValues]
asnIntNorm = [itm[1] for itm in asnValues]

aspIntAvg = [itm[0] for itm in aspValues]
aspIntNorm = [itm[1] for itm in aspValues]

N = len(asnIntAvg)
ind = np.arange(N)
width = 0.35

rects1 = axarr[0].bar(ind, asnIntAvg, width, color='r')
rects2 = axarr[0].bar(ind+width, asnIntNorm, width, color='y')


axarr[0].set_ylabel('Intensity Value', fontsize=40)
axarr[0].set_xlabel('Aspargine', fontsize=40)
axarr[0].set_xticks(ind+width)
axarr[0].set_xticklabels(asnAtoms,fontsize=20)

N = len(aspIntAvg)
ind = np.arange(N)

rects1 = axarr[1].bar(ind, aspIntAvg, width, color='r')
rects2 = axarr[1].bar(ind+width, aspIntNorm, width, color='y')

axarr[1].set_ylabel('Intensity Value', fontsize=40)
axarr[1].set_xlabel('Aspartic Acid', fontsize=40)
axarr[1].set_xticks(ind+width)
axarr[1].set_xticklabels(aspAtoms,fontsize=20)



fig.set_figheight(20)
fig.set_figwidth(40)
plt.savefig("AspVsAsn.png")



# plt.show()
