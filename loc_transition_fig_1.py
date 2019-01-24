import networkx as nx
import numpy as np
from numpy import linalg as LA
import random
import scipy as sp
import scipy.sparse as sparse
import matplotlib.pylab as plt

# contact information regarding codes
# "angeliya c.u" <008.angeliya@gmail.com>, "priodyuti pradhan" <priodyutipradhan@gmail.com>
# Complex Systems Lab, Indian Institute of Technology Indore 

#---------modelling network------------ 
n = 100 # Number of nodes
c = 3    # average degree of the random subgraph
deg = []
hub = {}
lst = []

p = c/float(n-1)   
print 'connection probability: ', p

flag = True
while flag:
  gnp = nx.fast_gnp_random_graph(n,p)
  if nx.is_connected(gnp):
     flag = False 
  
print 'Number of nodes: ',n
print 'Number of edges: ', len(gnp.edges())
print 'Check for connectedness: ',nx.is_connected(gnp)
A = nx.to_scipy_sparse_matrix(gnp,nodelist=None,dtype=np.float32,format='lil')
d = A.sum(axis=1)
avg_deg = sum(d)/n
print 'avg deg: ' ,avg_deg[0,0] 

#-------for loop for iterating d---------

fd = open('hub-deg_vs_ipr.txt','w')

for d in range(10,200,1):
 #fd1 = open(str(d)+'.txt','w') #To store the network in each step
 A[:,(n-1)] = 0
 A[(n-1),:] = 0
 deg.append(d)
 p_d = d/float(n-1)
 for i in range(n): 
    r = random.uniform(0,1)
    if r<p_d: 
       A[i,n-1] = 1
       A[n-1,i] = 1
 A[n-1,n-1] = 0
 hub_node = sum(A[:,(n-1)])
 #print 'size of the hub node: ', hub_node[0,0]

#---------IPR------------------------
 e_val, e_vec = sparse.linalg.eigsh(A,k=1)

 IPR = 0.0;
 for i in range(0,n):
   IPR = IPR + pow(e_vec[i],4)

 del e_vec
 fd.write("%d %f\n" % (int(d),float(IPR)))
 hub[d] = IPR
 lst.append(IPR)
 
 # To store the network in each step
 '''for p in range(0,n-1,1):
   for q in range(p+1,n,1):
      if A[p,q]==1:
         fd1.write(str(p+1)+' '+str(q+1)+'\n')
 fd1.close()'''              
  
fd.close() 

plt.plot(deg,lst,'ro--',markersize=2)
plt.grid()
plt.ylabel('IPR')
plt.xlabel('d')
plt.title('model network: localisation transition')
plt.show()




