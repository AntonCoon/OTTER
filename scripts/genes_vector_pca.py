#!/usr/bin/env python
# coding: utf-8

import sys
import os
import subprocess as sp
from pathlib import Path

import importlib
import scripts.settings as settings
importlib.reload(settings)


patient_folder = sys.argv[1]
if patient_folder[-1] != '/':
    patient_folder += '/'


vcfs = []
for f in os.listdir(patient_folder):
    if f.endswith('.vcf'):
        vcfs.append(f)

path_with_genes = settings.TMP_FOLDER + 'genes/'

run_anno = True
if run_anno:
    os.makedirs(path_with_genes, exist_ok=True)
    sp.run(f"rm -f {path_with_genes}/*", shell=True, check=True)

    for vcf in vcfs:
        fp = patient_folder + vcf
        out = settings.TMP_FOLDER + 'anno.vcf'
        cmd = f"java -Xmx8g -jar bin/snpEff/snpEff.jar GRCh38.105 -onlyProtein -no-intergenic -no-downstream -no-intron -no-upstream {fp} > {out}"
        sp.run(cmd, shell=True, check=True)
        sp.run(['mv', 'snpEff_genes.txt', path_with_genes + Path(vcf).stem + '_snpEff_genes.txt'], check=True)

# now goes the notebook.....

# In[84]:

###############################################3
###############################################3

import pandas
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from dataclasses import dataclass, field
from pathlib import Path
import numpy
import seaborn
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.patches as mpatches

seaborn.set_style("whitegrid")


# ## Parsing and preparation

# In[85]:


@dataclass(eq=True, order=True, frozen=True)
class Point:
    patient_id: str
    replica: str
    timepoint: str

@dataclass(eq=True, order=True, frozen=True)
class PointJoined:
    patient_id: str
    timepoint: str

# In[86]:



# In[87]:


data_tables = dict()
for genes_txt in os.listdir(path_with_genes):
    f = genes_txt.split('_')
    ts = f[1]
    rep = f[2][3:]

    point = Point(
        patient_id="0",
        replica=rep,
        timepoint=ts
    )
    data_tables[point] = pandas.read_csv(path_with_genes + genes_txt, sep="\t", header=1)


joined_data_tables = dict()

for patient_id in ['0']:
    for timepoint in ["1", "2", "3", "4", "5"]:
        joined_point = PointJoined(patient_id=patient_id, timepoint=timepoint)
        data_frame_rep1 = data_tables[Point(patient_id=patient_id, timepoint=timepoint, replica="1")]
        data_frame_rep2 = data_tables[Point(patient_id=patient_id, timepoint=timepoint, replica="2")]
        joinsed_data_frame = pandas.concat([data_frame_rep1, data_frame_rep2])

        joined_data_tables[joined_point] = joinsed_data_frame



# In[89]:


genes_set = set()
for point, data_frame in joined_data_tables.items():
    genes_set = genes_set.union(set(data_frame["GeneId"]))


# In[90]:


vector_size = len(genes_set)
genes_to_idx = {g_name: idx for idx, g_name in enumerate(genes_set)}



# In[92]:


points = list()
ordered_patients_points = []
for point, data_frame in joined_data_tables.items():

    # account for possible missing values
    for c in ['variants_impact_HIGH', 'variants_impact_MODERATE']:
        if c not in data_frame.columns:
            data_frame[c] = 0

    ordered_patients_points.append(point)
    gene_vector = numpy.zeros(vector_size)
    for gene_name, idx in genes_to_idx.items():
        gene_presence = numpy.sum(
            data_frame[
                data_frame["GeneId"] == gene_name
            ][
                ["variants_impact_HIGH", "variants_impact_MODERATE"]
            ].to_numpy()
        )
        if gene_presence:
            gene_vector[idx] = 1

    points.append(gene_vector)

# In[93]:


points = numpy.array(points)


# In[94]:



# ## Embeddings

# In[95]:


epochs_number = 4000


# In[96]:


class Autoencoder(nn.Module):
    def __init__(self, input_dim, hidden_dim, embedding_dim):
        super(Autoencoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, embedding_dim)
        )
        self.decoder = nn.Sequential(
            nn.Linear(embedding_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, input_dim),
            nn.Sigmoid()
        )

    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded


# In[99]:


masking_probabilities = numpy.random.normal(0.2, 0.01, epochs_number)
masking_probabilities[masking_probabilities < 0] = 0
masking_probabilities[masking_probabilities > 1] = 1


# In[100]:


data = torch.from_numpy(points).float()


# In[101]:


input_dim = data.shape[1]
hidden_dims = 128
embedding_dim = 8


# In[102]:


model = Autoencoder(input_dim, hidden_dims, embedding_dim)
criterion = nn.MSELoss()
# model.load_state_dict(torch.load("models/autoencoder_unmasker.pth"))
loss_values = []
optimizer = optim.Adam(model.parameters(), lr=0.003)
for epoch in range(epochs_number):
    # if (epoch + 1) % 100 == 0:
    #     print(f"Epoch [{epoch + 1}/{epochs_number}], Loss: {loss.item():.4f}")

    mask_prob = masking_probabilities[epoch]

    mask = torch.rand_like(data) > mask_prob
    masked_data = data * mask

    outputs = model(masked_data)
    loss = criterion(outputs, data)

    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    loss_values.append(loss.item())


# In[103]:


embeddings = model.encoder(data).detach().numpy()


# In[104]:


pca = PCA(n_components=2)
pca.fit(embeddings)


# In[105]:


components = pca.fit_transform(embeddings)

points_2d_dict = {
    "pca0": components.T[0],
    "pca1": components.T[1],
    "patient": [patient.patient_id for patient in ordered_patients_points],
    "timepoint": [patient.timepoint for patient in ordered_patients_points]
}

data_2d = pandas.DataFrame.from_dict(points_2d_dict)


# In[106]:
# fig, ax = plt.subplots(figsize=(8, 7))
# fig = seaborn.scatterplot(data=data_2d, x="pca0", y="pca1", hue="patient", style="patient")
# # ax.scatter(data_2d['pca0'], data_2d['pca1'])
#
# for line in range(0, data_2d.shape[0]):
#      plt.text(
#          data_2d.pca0[line]+0.2,
#          data_2d.pca1[line],
#          data_2d.timepoint[line],
#          horizontalalignment='left',
#          size='medium',
#          color='black',
#          weight='semibold'
#      )
#


# ## Flow calculation

# In[59]:


mean_x_for_patient_for_timepoint = {}
mean_y_for_patient_for_timepoint = {}
for point in ordered_patients_points:
    patient_id, timepoint = point.patient_id, point.timepoint
    subdata = data_2d[numpy.logical_and(data_2d["patient"] == patient_id, data_2d["timepoint"] == timepoint)]
    x_mean = numpy.mean(subdata["pca0"])
    y_mean = numpy.mean(subdata["pca1"])
    mean_x_for_patient_for_timepoint[(patient_id, timepoint)] = x_mean
    mean_y_for_patient_for_timepoint[(patient_id, timepoint)] = y_mean


# In[60]:


coord_x = []
coord_y = []
velocity_x = []
velocity_y = []
for idx, point in enumerate(ordered_patients_points):
    x, y = components[idx]
    v_x, v_y = 0, 0
    if (point.patient_id, str(int(point.timepoint) + 1)) in mean_x_for_patient_for_timepoint:
        v_x = mean_x_for_patient_for_timepoint[(point.patient_id, str(int(point.timepoint) + 1))] - x
        v_y = mean_y_for_patient_for_timepoint[(point.patient_id, str(int(point.timepoint) + 1))] - y
    coord_x.append(x)
    coord_y.append(y)
    velocity_x.append(v_x)
    velocity_y.append(v_y)



# In[61]:


x = numpy.array(coord_x)
y = numpy.array(coord_y)
dx = numpy.array(velocity_x)
dy = numpy.array(velocity_y)

plt.rcParams['image.cmap'] = 'Paired'

fig, ax = plt.subplots(figsize =(9, 9))

ax.quiver(
    x,
    y,
    dx,
    dy,
    [int(patient.patient_id) for patient in ordered_patients_points],
    units='width'
)
plt.close()


# ## Points bootstrap

# In[77]:


n = 155
bootstrapped_data = []
for mask_prob in masking_probabilities[:n]:
    mask = torch.rand_like(data) > mask_prob * 3
    bootstrapped_data.append(data * mask)


# In[78]:


bootstrapped_embeddings = model.encoder(torch.cat(bootstrapped_data, 0)).detach().numpy()


# In[79]:


pca = PCA(n_components=2)
pca.fit(bootstrapped_embeddings)



# In[80]:


bootstrapped_components = pca.fit_transform(bootstrapped_embeddings)

bootstrapped_points_2d_dict = {
    "pca0": bootstrapped_components.T[0],
    "pca1": bootstrapped_components.T[1],
    "patient": [patient.patient_id for patient in ordered_patients_points] * n,
    "timepoint": [patient.timepoint for patient in ordered_patients_points] * n
}

bootstrapped_data_2d = pandas.DataFrame.from_dict(bootstrapped_points_2d_dict)


# In[81]:


# for line in range(0, bootstrapped_data_2d.shape[0]):
#      print(
#          bootstrapped_data_2d.pca0[line]+0.2,
#          bootstrapped_data_2d.pca1[line],
#          bootstrapped_data_2d.timepoint[line]
#      )


# In[ ]:





# In[68]:


f, ax = plt.subplots(figsize=(8, 7))
fig = seaborn.scatterplot(data=bootstrapped_data_2d, x="pca0", y="pca1", hue="patient", style="patient")

for line in range(0, bootstrapped_data_2d.shape[0]):
     plt.text(
         bootstrapped_data_2d.pca0[line]+0.2,
         bootstrapped_data_2d.pca1[line],
         bootstrapped_data_2d.timepoint[line],
         horizontalalignment='left',
         size='medium',
         color='black',
         weight='semibold'
     )

os.makedirs(patient_folder + 'pipeline-out', exist_ok=True)
plt.savefig(patient_folder + 'pipeline-out/pca.pdf')
plt.close()

# skip following cells..........
sys.exit()

# ## Bootstraped flow

# In[82]:


def draw_flow():
    mean_x_for_patient_for_timepoint = {}
    mean_y_for_patient_for_timepoint = {}
    for point in ordered_patients_points:
        patient_id, timepoint = point.patient_id, point.timepoint
        subdata = bootstrapped_data_2d[numpy.logical_and(
            bootstrapped_data_2d["patient"] == patient_id,
            bootstrapped_data_2d["timepoint"] == timepoint
        )]
        x_mean = numpy.mean(subdata["pca0"])
        y_mean = numpy.mean(subdata["pca1"])
        mean_x_for_patient_for_timepoint[(patient_id, timepoint)] = x_mean
        mean_y_for_patient_for_timepoint[(patient_id, timepoint)] = y_mean

    coord_x = []
    coord_y = []
    velocity_x = []
    velocity_y = []
    for idx, point in enumerate(ordered_patients_points * n):
        x, y = bootstrapped_components[idx]
        v_x, v_y = 0, 0
        if (point.patient_id, str(int(point.timepoint) + 1)) in mean_x_for_patient_for_timepoint:
            v_x = mean_x_for_patient_for_timepoint[(point.patient_id, str(int(point.timepoint) + 1))] - x
            v_y = mean_y_for_patient_for_timepoint[(point.patient_id, str(int(point.timepoint) + 1))] - y
        coord_x.append(x)
        coord_y.append(y)
        velocity_x.append(v_x)
        velocity_y.append(v_y)

    x = numpy.array(coord_x)
    y = numpy.array(coord_y)
    dx = numpy.array(velocity_x)
    dy = numpy.array(velocity_y)

    plt.rcParams['image.cmap'] = 'tab10'

    fig, ax = plt.subplots(figsize =(10, 10))

    qv1 = ax.quiver(
        x,
        y,
        dx * 2,
        dy * 2,
        [int(patient.patient_id) for patient in ordered_patients_points] * n,
        units='width'
    )
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    handles = []
    for idx, manual_lable in enumerate(["1", "2", "3", "4"]):
        handles.append(mpatches.Patch(color=matplotlib.colormaps['tab10'](idx * 3), label=manual_lable))

    plt.legend(handles=handles, loc="upper center", ncol=4, title="Patients")

    plt.savefig("../results/patients_dynamic.pdf", bbox_inches = 'tight')
    plt.show()


# In[83]:


draw_flow()


# In[ ]:




