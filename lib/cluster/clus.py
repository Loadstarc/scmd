#!usr/bin/python3
"""
Agglomerative clustering with and without structure
===================================================

This example shows the effect of imposing a connectivity graph to capture
local structure in the data. The graph is simply the graph of 20 nearest
neighbors.

Two consequences of imposing a connectivity can be seen. First clustering
with a connectivity matrix is much faster.

Second, when using a connectivity matrix, average and complete linkage are
unstable and tend to create a few clusters that grow very quickly. Indeed,
average and complete linkage fight this percolation behavior by considering all
the distances between two clusters when merging them. The connectivity
graph breaks this mechanism. This effect is more pronounced for very
sparse graphs (try decreasing the number of neighbors in
kneighbors_graph) and with complete linkage. In particular, having a very
small number of neighbors in the graph, imposes a geometry that is
close to that of single linkage, which is well known to have this
percolation instability.
"""
# Authors: Gael Varoquaux, Nelle Varoquaux
# License: BSD 3 clause

import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph


# Generate sample data
# Create a graph capturing local connectivity. Larger number of neighbors
# will give more homogeneous clusters to the cost of computation
# time. A very large number of neighbors gives more evenly distributed
# cluster sizes, but may not impose the local manifold structure of
# the data

def cluster(nclus):
    np.random.seed(0)
    n_clusters_list  = [nclus]
    Y=np.loadtxt('cluster/tmpco')
    X = Y

    #print colors
    knn_graph = kneighbors_graph(X, 30, include_self=False)
    connectivity = None
    linkage = 'ward'
    for n_clusters in n_clusters_list :

            n_colors = n_clusters
            cm = plt.get_cmap('gist_rainbow', n_clusters)
            colors = []
            for i in range(n_colors) :
                color = cm(1.*i/n_colors)
                colors.append(color)
     
            model = AgglomerativeClustering(linkage=linkage,connectivity=connectivity,n_clusters=n_clusters)
            t0 = time.time()
            model.fit(X)
            elapsed_time = time.time() - t0
            
            fig, ax = plt.subplots()
            for k, col in zip(range(n_clusters), colors):
                my_members = model.labels_ == k
                ax.plot(X[my_members, 0], X[my_members, 1], 'w',
                     markerfacecolor=col, marker='o')

            plt.title('linkage=%s (time %.2fs)' % (linkage, elapsed_time),
                      fontdict=dict(verticalalignment='top'))

    #plt.savefig("test.png")
    np.savetxt('cluster/cluster_label.dat', model.labels_)
    #np.savetxt('colors')
