#!/usr/bin/env python 

import sys
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

matplotlib.use('Agg')

if len(sys.argv) < 3:
    print("Usage: infile outprefix")
    sys.exit(0)

file = sys.argv[1]
outprefix = sys.argv[2]

print("parsing "+file)

df = pd.read_csv(file, sep='\t')
print(df.head())
df.set_index(df.columns[0],inplace=True)
#sys.exit()
print("Read matrix: {} x {}".format(df.shape[0], df.shape[1]))
print(df)
#df.set_index('-')
print("Correlating..")
cormat = df.corr()
cormat = cormat.fillna(0)

print(cormat)

print("Performing decomposition.")
pca = PCA()
pca.fit(cormat)
components = pca.components_
expvar2 = pca.explained_variance_
explVar = pca.explained_variance_ratio_

print("{} - {}".format(len(expvar2),len( explVar)))

fhvar = open('explainedVariance.txt','w')
fhvar.write("PC\tExplainedVariance\tProportionalExplainedVariance\tCumulativeExplainedVariance\n")
cumulative = 0
print(expvar2)
for i in range(0,len(expvar2)):
	pc = i + 1
	cumulative += explVar[i]
	fhvar.write(f"PC{pc}\t{expvar2[i]}\t{explVar[i]}\t{cumulative}\n")
    # fhvar.write('PC' + str(pc))
fhvar.close()  

#                       print(explVar)
pcadf = pd.DataFrame(components)
pcadf.columns = cormat.columns
#                       pcadf.index.name = 'Component'
pcadf = pcadf.transpose()
colnames = []
for comp in range(len(components)):
    colnames.append("PC"+str(comp + 1))
pcadf.columns = colnames
#                       print(pcadf)

#pcadfmelt = pcadf.melt()
#                       print(pcadfmelt)

#pcadf.set_index('sample', inplace=True)
pcadf.index.name = "Sample"
print(pcadf.head())
pcadf = pcadf[pcadf.columns[:150]]
#print( pcadf.columns)
pcadf.to_csv(outprefix+"_PCs.txt", sep='\t')

fig, ax = plt.subplots()
sns.scatterplot(data=pcadf, x="PC1", y="PC2")
fig.savefig(outprefix+"_PC1and2.png")
plt.close(fig)

sns.scatterplot(data=pcadf, x="PC1", y="PC4")
fig.savefig(outprefix+"_PC1and4.png")
plt.close(fig)
