
# coding: utf-8

# **Anexo IV. Analisis anatómico-funcional**
# 
# **José Pedro Manzano**
# 
# **Lenguaje: Python**

# En este anexo se sigue el procedimiento necesario para, a partir de las matrices de conectividad obtenidas de cada dataset, construir una red binaria de la que extraer las principales propiedades, tanto globales como locales.
# 
# Para no hacer el anexo demasiado extenso, las principales conclusiones extraídas de los resultados pueden encontrarse en la memoria del trabajo.

# In[88]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import community as cm #python-louvain
import statistics as st
import operator as op
import itertools
import os
import sys
from sklearn import preprocessing


# **Importación de los datos**
# 
# En esta primera parte, se cargarán las matrices de conectividad individuales de cada sujeto con el objetivo de obtener un promedio de 10 sujetos de cada modalidad. También mantendremos la matriz de un sujeto de cada modalidad para realizar algunas comparaciones:
# 
# - **Redes Anatómicos (DTI)** de *sujetos sanos (healthy patients), con Alzheimer* y *Autismo*.
# - **Redes Funcionales (fMRI)** de *sujetos sanos (healthy patients), con ADHD' y *Autismo*.

# In[89]:


# Define the pathfiles
userpath = "/Users/hose/Desktop/TFM_TECI/"   ##### Define your work directory here and keep the tree structure

healthy_dti_path = userpath + "/Datos/Seleccionados/dti_sano/datos/"
apoe4_dti_path = userpath + "/Datos/Seleccionados/dti_alzheimer/datos/"
autism_dti_path = userpath + "/Datos/Seleccionados/dti_autism/datos/"
autism_fmri_path = userpath + "/Datos/Seleccionados/fmri_autismo/datos/"
adhd_fmri_path = userpath + "/Datos/Seleccionados/fmri_adhd/datos/"
healthy_fmri_path = userpath + "/Datos/Seleccionados/fmri_sanos/datos/"

pathfiles = [healthy_dti_path, apoe4_dti_path, autism_dti_path, autism_fmri_path, adhd_fmri_path, healthy_fmri_path]

# NAME OF NODES
hd_healthy_dti = pd.read_csv(userpath + "/Datos/Seleccionados/dti_sano/nodos/1013090_DTI_region_names_full_file.txt", header=None)
hd_apoe4_dti = pd.read_csv(userpath + "/Datos/Seleccionados/dti_alzheimer/nodos/APOE-4_2_region_names_full_file.txt", header=None)
hd_apoe4_dti = hd_apoe4_dti[0:110]
hd_autism_dti = pd.read_csv(userpath + "/Datos/Seleccionados/dti_autism/nodos/ASD47B_DTI_region_names_full_file.txt", header=None)
hd_autism_fmri = pd.read_csv(userpath + "/Datos/Seleccionados/fmri_autismo/nodos/ASD83B_rsfMRI_region_names_full_file.txt", header=None)
hd_adhd_fmri = pd.read_csv(userpath + "/Datos/Seleccionados/fmri_adhd/nodos/KKI_1018959_region_names_full_file.txt", header=None)
hd_healthy_fmri = pd.read_csv(userpath + "/Datos/Seleccionados/fmri_sanos/nodos/Baltimore_5560_region_names_full_file.txt", header=None)

node_names = [hd_healthy_dti, hd_apoe4_dti, hd_autism_dti, hd_autism_fmri, hd_adhd_fmri, hd_healthy_fmri]

## FOR THE SINGLE-SUBJECT LEVEL
individuals = {}
ind_names = ["healthy_dti_ind", "apoe4_dti_ind", "autism_dti_ind", "autism_fmri_ind", "adhd_fmri_ind", "healthy_fmri_ind"]
             
## FOR THE GROUP-LEVEL
group_avg = {}
groupnames = ["healthy_dti_group", "apoe4_dti_group", "autism_dti_group", "autism_fmri_group", "adhd_fmri_group", "healthy_fmri_group"]


# In[90]:


# Load connectivity matrixes
n = 0
for i in pathfiles:
    files = os.listdir(i)
    aux = pd.read_csv(i+files[0], header=None, delim_whitespace=True)
    k = 0
    for j in files:
        if k<1:
            print("Individual subject of", ind_names[n], "is", j)
            df = pd.read_csv(i+j, header=None, delim_whitespace=True)
            np.fill_diagonal(df.values,1)   # make sure diagonal = 0
            normalized_df = preprocessing.normalize(df, norm='max')  # normalize
            #normalized_df.columns = node_names[n]
            newname = ind_names[n]
            individuals[newname] = normalized_df
        df = pd.read_csv(i+j, header=None, delim_whitespace=True)
        np.fill_diagonal(df.values,1)
        aux = aux + df
        np.fill_diagonal(aux.values,1)   # make sure diagonal = 0
        k = k + 1
    
    avg_df = aux/len(files)   # average matrix of groups subjects
    np.fill_diagonal(normalized_df,1)
    normalized_df = preprocessing.normalize(avg_df, norm='max')  # normalize 
    #normalized_df.columns = node_names[n]
    newname = groupnames[n]
    group_avg[newname] = normalized_df
    n = n + 1  #update individual counter


# Puede ser interesante representar las matrices de conectividad (normalizadas) para tener una primera idea del aspecto de nuestros datos.

# In[91]:


#BY SUBJECTS
get_ipython().run_line_magic('matplotlib', 'inline')
fig = plt.figure(figsize=(12, 12))
plt.subplot(331)
np.fill_diagonal(individuals[ind_names[0]],0)
plt.imshow(individuals[ind_names[0]], cmap='hot', interpolation='nearest')
plt.title(ind_names[0], fontsize=20)
plt.subplot(332)
np.fill_diagonal(individuals[ind_names[1]],0)
plt.imshow(individuals[ind_names[1]], cmap='hot', interpolation='nearest')
plt.title(ind_names[1], fontsize=20)
plt.subplot(333)
plt.imshow(individuals[ind_names[2]], cmap='hot', interpolation='nearest')
plt.title(ind_names[2], fontsize=20)
plt.subplot(334)
plt.imshow(individuals[ind_names[3]], cmap='hot', interpolation='nearest')
plt.title(ind_names[3], fontsize=20)
plt.subplot(335)
plt.imshow(individuals[ind_names[4]], cmap='hot', interpolation='nearest')
plt.title(ind_names[4], fontsize=20)
plt.subplot(336)
plt.imshow(individuals[ind_names[5]], cmap='hot', interpolation='nearest')
plt.title(ind_names[5], fontsize=20)
plt.savefig('/Users/hose/Desktop/TFM_TECI/simulated_data/healthy_ill_ind')


# In[92]:


#BY AVERAGED GROUPS
fig2 = plt.figure(figsize=(12, 12))
plt.subplot(331)
plt.imshow(group_avg[groupnames[0]], cmap='hot', interpolation='nearest')
plt.title(groupnames[0], fontsize=20)
plt.subplot(332)
np.fill_diagonal(group_avg[groupnames[1]],0)
plt.imshow(group_avg[groupnames[1]], cmap='hot', interpolation='nearest')
plt.title(groupnames[1], fontsize=20)
plt.subplot(333)
plt.imshow(group_avg[groupnames[2]], cmap='hot', interpolation='nearest')
plt.title(groupnames[2], fontsize=20)
plt.subplot(334)
plt.imshow(group_avg[groupnames[3]], cmap='hot', interpolation='nearest')
plt.title(groupnames[3], fontsize=20)
plt.subplot(335)
plt.imshow(group_avg[groupnames[4]], cmap='hot', interpolation='nearest')
plt.title(groupnames[4], fontsize=20)
plt.subplot(336)
plt.imshow(group_avg[groupnames[5]], cmap='hot', interpolation='nearest')
plt.title(groupnames[5], fontsize=20)
plt.savefig('/Users/hose/Desktop/TFM_TECI/simulated_data/healthy_ill_group')


# Parece evidente que las matrices de conectividad promedio son capaces de filtrar gran parte de la variabilidad individual en las matrices funcionales, lo que revela de forma más evidente la estructura subyacente de la red. 
# 
# En cuanto a las redes anatómicas parece no existir la misma variabilidad que en las funcionales, ya que el promedio grupal y las individuales se mantienen bastante similares.

# **CONSTRUCCIÓN DEL GRAFO**
# 
# El objetivo final es comparar las propiedades de las diferentes redes y discutir si la teoría de grafos puede ser una herramienta útil o complementaria para detectar singularidades entre patologías.
# 
# Recuerda que las redes pueden contener valores negativos de correlación. Dado que no existe consenso sobre lo que implica una correlación negativa en DTI o en fMRI, se trabajará con los valores positivos en este caso.
# 

# In[93]:


np.any(individuals[ind_names[4]]<0)  # Any value below of 0? 


# In[94]:


individuals[ind_names[4]][individuals[ind_names[4]] < 0] = 0
np.any(individuals[ind_names[4]]<0)


# In[96]:


for i in range(0,6):
    individuals[ind_names[i]][individuals[ind_names[i]] < 0] = 0 
    group_avg[groupnames[i]][group_avg[groupnames[i]]< 0] = 0


# **DEFINE THE NULL-MODELS**
# 
# A continuación se definen como modelos nulos el modelo de Erdos-Renyi, una red en estrella y una red regular con propiedades similares a las de las redes construidas.

# In[97]:


n=110    # Erdos-Rengi nodes
m=1200   # Erdos-Rengi links
Grafos_ind = []
Grafos_group = []

# Creation of random graphs to compare with the brain graphs
Grandom=nx.gnm_random_graph(n,m,seed=123456789) #generador de grafo er, semilla fija
Gregular=nx.grid_2d_graph(10,11) #generador de grafo regular
Gstar=nx.star_graph(n-1) #generador de grafo en estrella

Grafos_ind.append(Grandom)
Grafos_ind.append(Gregular)
Grafos_ind.append(Gstar)

Grafos_group.append(Grandom)
Grafos_group.append(Gregular)
Grafos_group.append(Gstar)


# Como se indicaba en la memoria, como cada matriz proviene de un dataset diferente, tienen un número de nodos y enlaces diferente. Se intentará solventar este problema lo máximo posible para cada modalidad, ajustando el umbral de modo que las redes de DTI tengan  un número de enlaces y nodos lo más similar posible. Ídem para las de fMRI.

# In[98]:


#which = lambda lst:list(np.where(lst)[0])  # like which of R language
def num_of_zeros(xx,thr):
    np.fill_diagonal(x,0)
    xx[xx > thr] = 1
    xx[xx < thr] = 0
    current_links = np.sum(xx)#.values.sum()
    return(current_links)   #links different of 0


# In[99]:


# To relabel the nodes of the graph
def dict_names(x):
    dic={}
    lista_nombres = x[0].tolist()
    lista_nodos = range(len(lista_nombres))
    dic[lista_nodos]=lista_nombres
    return(dic)


# Y obtenemos las matrices de conectividad binarias tanto paras el promedio grupal como para el individual en cada modalidad:

# In[101]:


## FOR THE SINGLE-SUBJECT LEVEL
thr_ind = []
mat_ind = []
n = 0

for subject_i in individuals:
    x = individuals[subject_i]
    thr = 0.1
    print("Next subject!!!",subject_i)
    flag_error = 0
    flag_if = 0
    flag_else = 0
    while (abs(m - num_of_zeros(x,thr))>150 and (flag_error<2) and (thr<0.7)):       # while links in x > m, increase the threshold -- thr<0.7 stopping condition
        flag_error = flag_if + flag_else  # to avoid starting to jump into the if and the else (infinite loop)
        if (m < num_of_zeros(x,thr)):
            thr = thr + 0.03
            flag_if = 1
        else:
            thr = thr - 0.03
            flag_else = 1
    x[x > thr] = 1
    x[x < thr] = 0
    matrizAdy = np.matrix(x)
    mat_ind.append(matrizAdy)
    thr_ind.append(thr)
    G_ind = nx.from_numpy_matrix(matrizAdy)
    G_ind.name = "Graph_" + subject_i
    nx.relabel_nodes(G_ind, dict_names(node_names[n]))
    Grafos_ind.append(G_ind)
    n = n + 1


# In[102]:


## FOR THE GROUP LEVEL
thr_group = []
mat_group = []
n = 0

for group_i in group_avg:
    x = group_avg[group_i]
    thr = 0.1
    print("Next group average!!!", group_i)
    flag_error = 0
    flag_if = 0
    flag_else = 0
    while (abs(m - num_of_zeros(x,thr))>150 and (flag_error<2) and (thr<0.7)):       # while links in x > m, increase the threshold -- thr<0.7 stopping condition
        flag_error = flag_if + flag_else  # to avoid starting to jump into the if and the else (infinite loop)
        if (m < num_of_zeros(x,thr)):
            thr = thr + 0.03
            flag_if = 1
        else:
            thr = thr - 0.03
            flag_else = 1
    x[x > thr] = 1
    x[x < thr] = 0
    
    matrizAdy = np.matrix(x)
    mat_group.append(matrizAdy)
    
    thr_group.append(thr)
    G_group = nx.from_numpy_matrix(matrizAdy)
    G_group.name = "Graph_" + group_i
    nx.relabel_nodes(G_group, dict_names(node_names[n]))
    Grafos_group.append(G_group)
    n = n + 1


# In[103]:


plt.imshow(group_avg[group_i])


# In[104]:


group_avg[group_i]


# **ANÁLSIS I I - CARACTERÍSTICAS GLOBALES DE LA RED**

# In[105]:


## INDIVIDUAL SUBJECTS TABLE
tabla1 = []
tabla2 = []
tabla3 = []
for i in range(3,len(Grafos_ind)):
    tabla1.append(len(nx.nodes(Grafos_ind[i])))
    tabla2.append(len(nx.edges(Grafos_ind[i])))
    tabla3.append(thr_ind[i-3])

columnas=['Healthy_DTI','Alzheimer_DTI','Autism_DTI','Autism_fMRI','ADHD_fMRI','Healthy_fMRI']
ind = {0:'Nodes',1:'Edges',2:'Threshold'}
tabla= pd.DataFrame([tabla1,tabla2,tabla3],columns=columnas)
tabla.rename(index=ind)


# In[106]:


## GROUP TABLE
tabla1 = []
tabla2 = []
tabla3 = []
for i in range(3,len(Grafos_group)):
    tabla1.append(len(nx.nodes(Grafos_group[i])))
    tabla2.append(len(nx.edges(Grafos_group[i])))
    tabla3.append(thr_group[i-3])
    
columnas=['Healthy_DTI','Alzheimer_DTI','Autism_DTI','Autism_fMRI','ADHD_fMRI','Healthy_fMRI']
ind = {0:'Nodes',1:'Edges',2:'Threshold'}
tabla= pd.DataFrame([tabla1,tabla2,tabla3],columns=columnas)
tabla.rename(index=ind)


# In[107]:


# CARACHTERISTICS OF INDIVIDUAL-SUBJECTS GRAPH
Gcluster_ind=[]
Gdegree_ind=[]
Gbet_ind=[]
Gclo_ind=[]
Geff_ind=[]
Gasor_ind=[]
Gseg_ind=[]
Gint_ind=[]
Gpath_ind=[]

for i in range(len(Grafos_ind)):
    aux_degree = 0
    for v in Grafos_ind[i].nodes():
        aux_degree= aux_degree + Grafos_ind[i].degree(v)
    Gcluster_ind.append(nx.average_clustering(Grafos_ind[i]))
    Gdegree_ind.append((aux_degree)/(len(Grafos_ind[i].nodes())))
    Gbet_ind.append(st.mean(list(nx.betweenness_centrality(Grafos_ind[i]).values())))
    Gclo_ind.append(st.mean(list(nx.closeness_centrality(Grafos_ind[i]).values())))
    Geff_ind.append(nx.global_efficiency(Grafos_ind[i]))
    Gpath_ind.append(nx.average_shortest_path_length(Grafos_ind[i]))
    Gasor_ind.append(nx.degree_assortativity_coefficient(Grafos_ind[i]))
    Gseg_ind.append(Gcluster_ind[i]/Gcluster_ind[0]) #se considera el grafo e-r
    Gint_ind.append(Gpath_ind[i]/Gpath_ind[0])#integración con la eficiencia(en lugar del path)


# In[108]:


# CARACHTERISTICS OF AVERAGE GROUPS GRAPH

Gcluster_group=[]
Gdegree_group=[]
Gbet_group=[]
Gclo_group=[]
Geff_group=[]
Gasor_group=[]
Gseg_group=[]
Gint_group=[]
Gpath_group=[]

for i in range(len(Grafos_group)):
    aux_degree = 0
    for v in Grafos_group[i].nodes():
        aux_degree= aux_degree + Grafos_group[i].degree(v)
    Gcluster_group.append(nx.average_clustering(Grafos_group[i]))
    Gdegree_group.append((aux_degree)/(len(Grafos_group[i].nodes())))
    Gbet_group.append(st.mean(list(nx.betweenness_centrality(Grafos_group[i]).values())))
    Gclo_group.append(st.mean(list(nx.closeness_centrality(Grafos_group[i]).values())))
    Geff_group.append(nx.global_efficiency(Grafos_group[i]))
    Gpath_group.append(nx.average_shortest_path_length(Grafos_group[i]))
    Gasor_group.append(nx.degree_assortativity_coefficient(Grafos_group[i]))
    Gseg_group.append(Gcluster_group[i]/Gcluster_group[0]) #se considera el grafo e-r
    Gint_group.append(Gpath_group[i]/Gpath_group[0])#integración con la eficiencia(en lugar del path)


# Podemos mostrar las propiedades extraídas en una tabla para las redes individuales:

# In[109]:


# SUMMARY TABLE OF INDIVIDUAL GRAPHS
columnas=['Erdos-Renyi','Regular-Graph','Star-Graph','Sano_DTI','Alzheimer_DTI','Autismo_DTI','Autismo_fMRI','ADHD_fMRI','Sano_fMRI']
ind = {0:'Clustering',1:'Grado Medio',2:'Betweenness Centrality',3:'Closeness Centrality',4:'Global Efficiency',5:'Average Shortest Path',6:'Asortatividad',7:'Segregación(SW)',8:'Integración(SW)'}
tabla= pd.DataFrame([Gcluster_ind,Gdegree_ind,Gbet_ind,Gclo_ind,Geff_ind,Gpath_ind,Gasor_ind,Gseg_ind,Gint_ind],columns=columnas)
tabla.rename(index=ind)


# Y para las redes grupales:

# In[110]:


# SUMMARY TABLE OF GROUP GRAPHS
columnas=['Erdos-Renyi','Regular-Graph','Star-Graph','Sano_DTI','Alzheimer_DTI','Autismo_DTI','Autismo_fMRI','ADHD_fMRI','Sano_fMRI']
ind = {0:'Clustering',1:'Grado Medio',2:'Betweenness Centrality',3:'Closeness Centrality',4:'Global Efficiency',5:'Average Shortest Path',6:'Asortatividad',7:'Segregación(SW)',8:'Integración(SW)'}
tabla= pd.DataFrame([Gcluster_group,Gdegree_group,Gbet_group,Gclo_group,Geff_group,Gpath_group,Gasor_group,Gseg_group,Gint_group],columns=columnas)
tabla.rename(index=ind)


# **ANÁLISIS II - PROPIEDADES LOCALES DE LA RED**
# 
# De manera similar, podemos extraer las mismas propiedades para cada nodo en lugar de para toda la red.

# In[111]:


## FOR INDIVIDUAL SUBJECTS
Local_prox_ind = []
Local_bet_ind = []
Local_clo_ind = []
Local_cluster_ind = []
Grado_node_ind = []
for i in range(0,len(Grafos_ind)):
    Local_prox_ind.append(nx.average_neighbor_degree(Grafos_ind[i]))
    Local_bet_ind.append(nx.betweenness_centrality(Grafos_ind[i]))
    Local_clo_ind.append(nx.closeness_centrality(Grafos_ind[i]))
    Local_cluster_ind.append(nx.clustering(Grafos_ind[i]))
    grado_ind={}
    for v in Grafos_ind[i].nodes():
            grado_ind[v]=Grafos_ind[i].degree(v)
    Grado_node_ind.append(grado_ind)
    
for i in range(0,len(Grafos_ind)):
    Local_prox_ind[i] = sorted(Local_prox_ind[i].items(),reverse=True,key=op.itemgetter(1))
    Local_bet_ind[i] = sorted(Local_bet_ind[i].items(),reverse=True,key=op.itemgetter(1))
    Local_clo_ind[i] = sorted(Local_clo_ind[i].items(),reverse=True,key=op.itemgetter(1))
    Local_cluster_ind[i] = sorted(Local_cluster_ind[i].items(),reverse=True,key=op.itemgetter(1))
    Grado_node_ind[i] = sorted(Grado_node_ind[i].items(),reverse=True,key=op.itemgetter(1))


# In[112]:


## FOR GROUP AVERAGE SUBJECTS
Local_prox_group = []
Local_bet_group = []
Local_clo_group = []
Local_cluster_group = []
Grado_node_group = []
for i in range(0,len(Grafos_group)):
    Local_prox_group.append(nx.average_neighbor_degree(Grafos_group[i]))
    Local_bet_group.append(nx.betweenness_centrality(Grafos_group[i]))
    Local_clo_group.append(nx.closeness_centrality(Grafos_group[i]))
    Local_cluster_group.append(nx.clustering(Grafos_group[i]))
    grado_group={}
    for v in Grafos_group[i].nodes():
            grado_group[v]=Grafos_group[i].degree(v)
    Grado_node_group.append(grado_group)
    
for i in range(0,len(Grafos_group)):
    Local_prox_group[i] = sorted(Local_prox_group[i].items(),reverse=True,key=op.itemgetter(1))
    Local_bet_group[i] = sorted(Local_bet_group[i].items(),reverse=True,key=op.itemgetter(1))
    Local_clo_group[i] = sorted(Local_clo_group[i].items(),reverse=True,key=op.itemgetter(1))
    Local_cluster_group[i] = sorted(Local_cluster_group[i].items(),reverse=True,key=op.itemgetter(1))
    Grado_node_group[i] = sorted(Grado_node_group[i].items(),reverse=True,key=op.itemgetter(1))


# In[113]:


## SUMMARY TABLE
def local_properties(graph_x, propertie_x):
    indc={0:'Erdos-Renyi',1:'Regular-Graph',2:'Star-Graph',3:'Healthy_DTI',4:'Alzheimer_DTI',5:'Autism_DTI',6:'Autism_fMRI',7:'ADHD_fMRI',8:'Healthy_fMRI'}
    lista=[]
    for i in range(0,len(graph_x)):
        lista.append(itertools.islice(propertie_x[i],0,4))
    tabla_x= pd.DataFrame(lista)
    tabla_x.rename(index=indc)
    return(tabla_x)


# En este caso, resumir y mostrar la información es algo más complejo (como se vio antes, son cientos de nodos). Por ello, mostraremos el valor de los 4 nodos más extremos para cada propiedad, en forma de pares de valor (nodo, valor) para tener una idea:
# 
# **SHORT PATH LENGHTS**
# 
# *INDIVIDUAL SUBJECTS*

# In[61]:


indc={0:'Erdos-Renyi',1:'Regular-Graph',2:'Star-Graph',3:'Healthy_DTI',4:'Alzheimer_DTI',5:'Autism_DTI',6:'Autism_fMRI',7:'ADHD_fMRI',8:'Healthy_fMRI'}
columns = {}


# In[62]:


local_properties(Grafos_ind, Local_prox_ind).rename(index=indc)


# *GROUP AVERAGE*

# In[ ]:


local_properties(Grafos_group, Local_prox_group).rename(index=indc)


# **BETWENNESS**

# *INDIVIDUAL SUBJECTS*

# In[64]:


local_properties(Grafos_ind, Local_bet_ind).rename(index=indc)


# *GROUP AVERAGE*

# In[65]:


local_properties(Grafos_group, Local_bet_group).rename(index=indc)


# **CLOSENESS**

# *INDIVIDUAL SUBJECTS*

# In[66]:


local_properties(Grafos_ind, Local_clo_ind).rename(index=indc)


# *GROUP AVERAGE*

# In[67]:


local_properties(Grafos_group, Local_clo_group).rename(index=indc)


# **CLUSTERING COEFFICIENT**

# *INDIVIDUAL SUBJECTS*

# In[68]:


local_properties(Grafos_ind, Local_cluster_ind).rename(index=indc)


# *GROUP AVERAGE*

# In[69]:


local_properties(Grafos_group, Local_cluster_group).rename(index=indc)


# **DEGREE NODE**

# *INDIVIDUAL SUBJECTS*

# In[70]:


local_properties(Grafos_ind, Grado_node_ind).rename(index=indc)


# *GROUP AVERAGE*

# In[71]:


local_properties(Grafos_group, Grado_node_group).rename(index=indc)

