# Before run this code please run roots_of_cubic_equations_fig_4.m 
# for a given number of nodes and number of edges

# For this code:
# Input: k, n1, n2
# k is the average degree of the random regular graph
# n1 is the size of wheel graph
# n2 is the size of random regular graph
# k, n1, n2 obtain from roots_of_cubic_equations_fig_4.m
# k can be k1 or k2 such that k1*n2 or k2*n2 should be even 

#Output:1) one file stored degree sequence 
#       2) one file stored principal eigenvector entries
#       3) one file stored the wheel-random-regular graph 


# contact information regarding codes
# "priodyuti pradhan" <priodyutipradhan@gmail.com>
# Complex Systems Lab, Indian Institute of Technology Indore 

'''
It simulates the optimized structure: two graph components (wheel and random regular) connected via a node.
'''
import networkx as nx                                              
import matplotlib.pyplot as plt 
import numpy as np
import scipy as sp
import codecs
import sys
from operator import itemgetter
from numpy import linalg as LA
from functions_base import store_undirected_network


def store_undirected_network(A, N, path):
  f_out = open(path,'wb')
  for i in range(N-1):
    for j in range(i+1,N,1):
      if A[i,j] == 1:
        #print i+1,j+1
        f_out.write(str(i+1)+' '+str(j+1)+'\n')
  f_out.close()

def get_ipr(A):
  e_val, e_vec = LA.eigh(A)
  #print e_val
  d = np.shape(e_vec)   
  n = d[0]   
  lam1 = e_val[n-1]  
  lam2 = e_val[n-2]   
  Ev1 = abs(e_vec[:,n-1])
  Ev2 = e_vec[:,n-2]
            
  IPR1 = 0.0;
  IPR2 = 0.0;
  for i in range(n):
    IPR1 = IPR1 + pow(Ev1[i],4);
    IPR2 = IPR2 + pow(Ev2[i],4); 
  path = 'wheel_RR_graph_pev_dloc_'+str(n)+'.txt'   

  f_out = open(path,'wb')
  for i in range(n): 
    f_out.write(str(i+1) + ' ' + str(Ev1[i])+'\n')  
  f_out.close()
  
  del e_val
  del e_vec
  del Ev1
  del Ev2
  return IPR1,IPR2,lam1,lam2


# k, n1, and n2 values are obtain from the roots_of_cubic_equations_fig_4.m

k = 45          
n1 = 1937       
n2 = 470         #  We will increase the n2 size by adding more number of nodes with same average degree k
                 # Which will automatically increase the size N 

print 's1: ',n1
Flag = False
while not Flag:
  G2 = nx.random_regular_graph(k,n2)
  Flag = nx.is_connected(G2)
  print Flag

G1 = nx.wheel_graph(n1)
deg_seq = nx.degree(G1)
max_deg = 0
node_index = 0
for i in list(nx.degree(G1).keys()):
  if deg_seq[i] > max_deg:
     max_deg = deg_seq[i];
     node_index = i

print 'max deg: ',max_deg
print 'Node label with max deg: ',node_index


N = n1 + n2 + 1
M = np.zeros((N,N), dtype=int)
M1 = np.zeros((n1,n1), dtype=int)
M2 = np.zeros((n2,n2), dtype=int)

Edge_list1 = nx.edges(G1)
for edge in Edge_list1:
  x = edge[0]
  y = edge[1]
  M1[x][y] = 1
  M1[y][x] = 1

Edge_list2 = nx.edges(G2)
for edge in Edge_list2:
  x = edge[0]
  y = edge[1]
  M2[x][y] = 1
  M2[y][x] = 1

#print 'Number of nodes in random regular graph: %d' %s1
Edges1 = nx.number_of_edges(G1)
avg_deg1 = Edges1*2/float(n1)
print 'Number of edges in wheel graph: ', Edges1
print 'Average degree of wheel graph : %f' %avg_deg1
#ipr1,ipr2,lam1,lam2 = get_ipr(M1)                
#print 's1: %d ipr1: %0.4f ipr2: %0.4f lam1: %0.16f lam2: %0.4f\n'%(s1,ipr1,ipr2,lam1,lam2)                                                                                                        


#print 'Number of nodes in wheel graph: %d' %s2
Edges2 = nx.number_of_edges(G2)
avg_deg2 = Edges2*2/float(n2)
print 'Number of edges in random regular graph: ', Edges2
print 'Average degree of random regular graph : %f' %avg_deg2
#ipr1,ipr2,lam1,lam2 = get_ipr(M2)                
#print 's2: %d ipr1: %0.4f ipr2: %0.4f lam1: %0.16f lam2: %0.4f\n'%(s2,ipr1,ipr2,lam1,lam2)                                                                                                        

for i in range(n1):
 for j in range(n1):
    M[i][j] = M1[i][j]
#print(len(M2))
#print(len(M))

M[n1-1][n1] = 1
M[n1][n1-1] = 1 
M[n1][n1+1] = 1
M[n1+1][n1] = 1

i = n1 + 1
j = n1 + 1
m = 0
n = 0

while i<N and m<n2:
   j = n1+1
   n = 0
   while j<N and n<n2:
     M[i][j] = M2[m][n]
     #print i,j,m,n
     j = j + 1
     n = n + 1
   i = i + 1
   m = m + 1

# For observing delocalization transition uncomment the below two lines
#M[0][6] = 0
#M[6][0] = 0

'''M[0][8] = 0
M[8][0] = 0
'''
'''for index in range(n1+4,n1+n2,1):
  if M[index][index+1]==0:  
      M[index][index+1] = 1
      M[index+1][index] = 1   
      break 
'''
gp = nx.from_numpy_matrix(M)
Edges = nx.number_of_edges(gp)
avg_deg = Edges*2/float(N)
print 'Number of nodes in wheel-random_regular_graph_combined: %d' %N
print 'Number of edges in wheel-random_regular_graph_combined: %d' %Edges
print 'Average degree of wheel-random_regular_graph_combined : %f' %avg_deg

ipr1,ipr2,lam1,lam2 = get_ipr(M)                
print 'N: %d \n ipr1: %0.16f \n ipr2: %0.16f \n lam1: %0.16f \n lam2: %0.16f\n'%(N,ipr1,ipr2,lam1,lam2)                                                                                                       
print(nx.is_connected(gp))
#print(gp.degree())

deg_vec = gp.degree()

dpath = 'wheel_RRd_deg_'+str(N)+'.txt'
f_deg = open(dpath,'wb')

for key, value in deg_vec.iteritems():
    f_deg.write(str(key+1) + ' ' + str(value)+'\n')
    #print key+1, value
f_deg.close()

path = 'wheel_RRd_graph_'+str(N)+'.txt'
store_undirected_network(M, N, path)


