import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
df= pd.read_excel('CVM-2020-NARMS-AnimalPathogenData.xlsx', nrows=14000)
print(len(df))



df = df.drop(df[df['Genus'] != 'E. coli'].index)
print(len(df))



print(df.iloc[400][:])
print('***************************')
print(df.iloc[900][:])


df_MIC = df.pivot(index = 'Sample ID', columns='Drug Name', values='MIC')
print(len(df['Sample ID'].unique()))
df_MIC.dropna(inplace=True)

df_year = df.pivot(index = 'Sample ID', columns='Drug Name', values='Year')
df_year.dropna(inplace=True)







print(df_year['Amikacin'].unique())


X = df_MIC.to_numpy()
model = TSNE(n_components=2, perplexity=5)
X_reduced = model.fit_transform(X)




scatter_x = X_reduced[:, 0]
scatter_y = X_reduced[:, 1]
group = np.array(df_year['Amikacin'])
cdict = {2017: 'orange', 2018: 'red', 2019: 'blue', 2020: 'green'}

fig, ax = plt.subplots()
for g in np.unique(group):
    ix = np.where(group == g)
    ax.scatter(scatter_x[ix], scatter_y[ix], c = cdict[g], label = int(g), s = 20)
ax.legend()
plt.show()