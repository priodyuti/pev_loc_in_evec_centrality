import networkx as nx
import numpy as np
from numpy import linalg as LA
import random
import scipy as sp
import scipy.sparse as sparse
import matplotlib.pylab as plt

# "angeliya c.u" <008.angeliya@gmail.com>, "priodyuti pradhan" <priodyutipradhan@gmail.com>
# Complex Systems Lab, Indian Institute of Technology Indore 
#-----------------------------------modelling network-------------------------------------------------
# vary n = 500, 1000, 10000, 50000, 100000, 150000 

n = 1000
c = 10
deg = []
hub = {}
lst = []
p = c/float(n-1)   
print('connection probability: ', p)

flag = True
while flag:
  gnp = nx.fast_gnp_random_graph(n,p)
  if nx.is_connected(gnp):
     flag = False 
  
print('Number of nodes: ',n)
print 'Number of edges: ', len(gnp.edges())
print('Check for connectedness: ',nx.is_connected(gnp))

A = nx.to_scipy_sparse_matrix(gnp,nodelist=None,dtype=np.float32,format='lil')
d = A.sum(axis=1)
avg_deg = sum(d)/n
print('avg deg: ' ,avg_deg[0,0]) 

#-------for loop for iterating d---------
path = 'c_'+str(n)+'_'+str(c)+'.txt'
fd = open(path,'w')
d = c*(c+1) + 30

A[:,(n-1)] = 0
A[(n-1),:] = 0

p_d = d/float(n-1)
#print(d)
for i in range(n): 
  r = random.uniform(0,1) 
  if r<p_d: 
     A[i,n-1] = 1
     A[n-1,i] = 1
A[n-1,n-1] = 0
hub_node = sum(A[:,(n-1)])
print('size of the hub node: ', hub_node[0,0])

#---------------IPR------------------------
e_val, e_vec = sparse.linalg.eigsh(A,k=1)


IPR = 0.0;
for i in range(0,n):
  IPR = IPR + pow(e_vec[i],4)

print("IPR = %f \n" % (float(IPR)))

e_vec = abs(e_vec)
Ev = sorted(e_vec,key = lambda x: x[0], reverse = True)
x=range(1,n+1,1)

plt.loglog(x,Ev, 'o--c', label = 'n = '+str(n))
plt.xlabel('i', fontsize = 15)
plt.ylabel(r'$ (v_1)_i $', fontsize = 15)
plt.legend(loc = 'best')
plt.grid(True)
plt.show()

index = 1
for i in range(n):
  fd.write(str(int(index))+ ' ' +str(float(Ev[i]))+ ' \n')
  index = index + 1
del e_vec
 
fd.close() 

