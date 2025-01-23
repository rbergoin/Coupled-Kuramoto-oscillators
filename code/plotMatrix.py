#!/usr/bin/env python

"""
Load and plot data saved

Developped by Raphael BERGOIN

Run : python3 plotMatrix.py
"""


#python3 -mpip install LIB_NAME --user
from math import *
import string
import numpy as np
from io import StringIO
from matplotlib.colors import ListedColormap
import codecs
import random
import warnings
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import copy
import cmath
import operator
from scipy import stats
from sklearn.metrics import mutual_info_score
import networkx as nx
import community as community_louvain
from sklearn.metrics import mutual_info_score
from scipy.stats import entropy


warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
    

def orderParameter(phi, m, N) :
    """
        Calculate the order parameter
        
        Parameters : 
        phi -- phase vector
        m -- order
        N -- number of neurons
    """
    
    R = 0.0
    for i in range(N) :
        R += cmath.exp(1.0j*m*phi[i])
        
    return abs((1.0/N)*R)
    
    
def localFields(phi, k, N) :
    """
        Calculate the local field of each neuron at time t
        Measure the degree of synchronization of neurons among others
        
        Parameters : 
        phi -- phase vector
        k -- weights matrix
        N -- number of neurons
    """
    
    z = np.zeros(N, dtype=complex)
    for i in range(N) :
        for j in range(N) :
            z[i] += k[i][j]*cmath.exp(1.0j*phi[j])
        
    return abs((1.0/N)*z)
    

def localOrderParameter(phi, N) :
    """
        Calculate the local order parameter between each pairs of neurons
        
        Parameters : 
        phi -- phase vector
        N -- number of neurons
    """
    
    localOrderParameters = np.zeros((N,N))
    
    for i in range(N) :
        for j in range(N) :
            localOrderParameters[i][j] = cos(phi[i]-phi[j])
        
    return localOrderParameters
    
    
def autocorrelation(phi, phi2, N) :
    """
        Autocorrelation function
        
        Parameters : 
        phi -- phase vector at time 1
        phi2 -- phase vector at time 2
        N -- number of neurons
    """
    
    C = 0.0
    for i in range(N) :
        C += cmath.exp(1.0j*phi2[i])*cmath.exp(-1.0j*phi[i])
    
    return np.dot(abs((1.0/N)*C), abs((1.0/N)*C))
    



""" Arguments """
    
if len(sys.argv) > 1 :
    if str(sys.argv[1]) == "1" :
        save = True
    else :
        save = False
else :
    save = False


"""Get data saved"""   

f = open("weights_matrices.txt", "r")
matrices = []
matrix = []
for x in f:
    if len(x.strip()) == 0 :
        matrices.append(matrix)
        matrix = []
    else :
        lst = x.split()
        matrix.append([float(i) for i in lst])

f.close()
matrices = np.array(matrices)


f = open("adjacency.txt", "r")
adjacency = []
for x in f:
    lst = x.split()
    adjacency.append([int(i) for i in lst])

f.close()
adjacency = np.array(adjacency)


f = open("phases.txt", "r")
phases = []
for x in f:
    lst = x.split()
    phases.append([float(i) for i in lst])

f.close()
phases = np.array(phases)


f = open("changeRates.txt", "r")
changeRates = []
for x in f:
    lst = x.split()
    changeRates = [float(i) for i in lst]

f.close()
changeRates = np.array(changeRates)


f = open("natural_frequencies.txt", "r")
naturalFrequencies = []
for x in f:
    lst = x.split()
    naturalFrequencies = [float(i) for i in lst]

f.close()
naturalFrequencies = np.array(naturalFrequencies)



orderParameter1 = []    #list of order parameters (m=1)
orderParameter2 = []    #list of order parameters (m=2)
orderParameter3 = []    #list of order parameters (m=3)
orderParameter4 = []    #list of order parameters (m=4)
localFields1 = []        #list of local fields through the time
autocorrelations = [0.0 for x in range(len(phases)//5)]   #list of autocorrelations


adimensional = True     #If adimensional time
dt = 0.1               #Time step
if adimensional :
    tau = 1.0
else : 
    tau = 0.01
dt = dt*tau
tfinal = (len(phases)-1.0)*dt
tpoints = np.arange(0.0, (tfinal+0.00001), dt)

ta = 0
#Calculate order parameter and autocorrelation from phase evolution
for t in range(0, len(phases)) :
    orderParameter1.append(orderParameter(phases[t], 1, len(phases[t])))   #calculate and register the order parameters (m=1)
    orderParameter2.append(orderParameter(phases[t], 2, len(phases[t])))   #calculate and register the order parameters (m=2)
    orderParameter3.append(orderParameter(phases[t], 3, len(phases[t])))   #calculate and register the order parameters (m=3)
    orderParameter4.append(orderParameter(phases[t], 4, len(phases[t])))   #calculate and register the order parameters (m=4)
    
    ta += 1
    if ta >= len(phases)/5 :
        autocorrelations.append(autocorrelation(phases[t-len(phases)//5], phases[t], len(phases[t])))  #calculate and register autocorrelations

#Moving average
orderParameter1 = np.convolve(orderParameter1, np.ones(300), 'valid') / 300
orderParameter2 = np.convolve(orderParameter2, np.ones(300), 'valid') / 300
orderParameter3 = np.convolve(orderParameter3, np.ones(300), 'valid') / 300
orderParameter4 = np.convolve(orderParameter4, np.ones(300), 'valid') / 300

#Calculate local fields evolution
for t in range(0, len(matrices)) :
    localFields1.append(localFields(phases[int(t*((len(phases)-1)/(len(matrices)-1)))], matrices[t], len(phases[t]))) #calculate and register the local fields of neurons
localFields1 = np.array(localFields1)



#Order according to phases
order = np.argsort(phases[-1,:])
    



newcmap = np.genfromtxt('./Matplotlib_colourmap_div.csv', delimiter=',')
cmap_div = ListedColormap(newcmap)





"""Plot data""" 

 

#Time development of the order parameters
plt.plot(tpoints[0:len(orderParameter1)], orderParameter1, label='order (m=1)') 
plt.plot(tpoints[0:len(orderParameter2)], orderParameter2, label='order (m=2)') 
plt.plot(tpoints[0:len(orderParameter3)], orderParameter3, label='order (m=3)') 
plt.plot(tpoints[0:len(orderParameter4)], orderParameter4, label='order (m=4)') 
#add orders if necessary
plt.gca().legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
if save : 
    plt.savefig('results/order_parameters.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Time development of the global order parameters', fontsize=25)
    plt.show() 
    
    
#Time development of the local fields
for i in range(len(localFields1[0])) :
    plt.plot(np.arange(0, tfinal+0.00001, dt*int((len(phases)-1)/(len(matrices)-1))), localFields1[:,i], label='Neuron %d' % i) 
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
#plt.gca().legend(fontsize=20)
if save : 
    plt.savefig('results/local_fields.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Time development of the local fields', fontsize=25)
    plt.show()
    
    
#Local order parameter at the end
localOrderParameters = localOrderParameter(phases[-1], len(phases[-1]))
plt.matshow(localOrderParameters, cmap=plt.cm.seismic, vmin=0, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
if save : 
    plt.savefig('results/local_order_parameters.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Local order parameters at the end', fontsize=25)
    plt.show() 


"""
#Adjacency matrix
#kij : j presynaptic to i postsynaptic neurons
plt.matshow(adjacency, cmap=plt.cm.gray_r, vmin=0, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
if save : 
    plt.savefig('results/adjacency_matrix.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Adjacency matrix', fontsize=25)
    plt.show() 
"""



#Weights matrix after updating (sorted) 
#kij : j presynaptic to i postsynaptic neurons
plt.matshow(matrices[0][order, :][:, order], cmap=cmap_div, vmin=-1, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
if save : 
    plt.savefig('results/weights_matrix_sorted_T0.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Weight matrix (sorted) T0', fontsize=25)
    plt.show() 
  
#Weights matrix after updating (sorted) 
#kij : j presynaptic to i postsynaptic neurons
plt.matshow(matrices[1][order, :][:, order], cmap=cmap_div, vmin=-1, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
if save : 
    plt.savefig('results/weights_matrix_sorted_T1.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Weight matrix (sorted) T1', fontsize=25)
    plt.show()

"""
#Weights matrix after updating (sorted) 
#kij : j presynaptic to i postsynaptic neurons
plt.matshow(matrices[2][order, :][:, order], cmap=cmap_div, vmin=-1, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
if save : 
    plt.savefig('results/weights_matrix_sorted_T2.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Weight matrix (sorted) T2', fontsize=25)
    plt.show()         
"""
        
#Weights matrix after updating (sorted) 
#kij : j presynaptic to i postsynaptic neurons
plt.matshow(matrices[-1][order, :][:, order], cmap=cmap_div, vmin=-1, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
if save : 
    plt.savefig('results/weights_matrix_sorted_T3.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Weight matrix (sorted) T3', fontsize=25)
    plt.show() 



#Distribution of the weights
plt.hist(matrices[0].flatten(), bins=100, histtype=u'step', density=True, label='before updating', log=True)
plt.hist(matrices[len(matrices)//2].flatten(), bins=100, histtype=u'step', density=True, color='green', label='in the middle', log=True)
plt.hist(matrices[-1].flatten(), bins=100, histtype=u'step', density=True, color='red', label='at the end', log=True)
plt.gca().set_xlabel('Weights', fontsize=20)
plt.gca().set_ylabel('Population', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.gca().legend(fontsize=20)
if save : 
    plt.savefig('results/distribution_weights.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    #plt.gca().set_title('Distribution of the weights', fontsize=25)
    plt.gca().set_title('Distribution of the weights (log scale)', fontsize=25)
    plt.show() 
    

#Evolution of the change rate of weights
plt.plot(tpoints[1:len(tpoints)], changeRates) 
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
if save : 
    plt.savefig('results/rateChange_weights.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Average evolution of the change rate of weights', fontsize=25)
    plt.show() 
    
    
"""
#Evolution  phases
#for k in range(len(phases[0])) :
    #plt.plot(tpoints, phases[:,k], label='Phase neuron %d' % k) 
plt.plot(tpoints, phases[:,0], label='Phase neuron 0') 
plt.plot(tpoints, phases[:,1], label='Phase neuron 1') 
#plt.plot(tpoints, np.subtract(phases[:,0], phases[:,1]), label='Difference of phases neuron 0 and 1') 
plt.gca().legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
if save : 
    plt.savefig('results/phases_evolution.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Evolution of  phases', fontsize=25)
    plt.show()   
"""    
    
    
#Strength of neurons
threshold = 0.0     #If necessary
a = np.absolute(matrices[-1])
#a = np.absolute(matrices[-1][order, :][:, order])
plt.plot(np.sum(np.where(a>threshold, a, 0), axis=1), '_', label='Incoming weights')
plt.plot(np.sum(np.where(a>threshold, a, 0), axis=0), '_', label='Outcoming weights', color='red')
plt.plot(np.sum(np.where(a>threshold, a, 0), axis=0)+np.sum(np.where(a>threshold, a, 0), axis=1), '_', label='Both weights', color='green')
plt.gca().set_ylim([0, len(matrices[-1])*2]) 
plt.gca().legend(fontsize=20)
plt.gca().set_ylabel('Strength', fontsize=20)
plt.gca().set_xlabel('Neurons', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/strength_neurons.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Strength of neurons', fontsize=25)
    plt.show() 
    
    
#Distribution strength
plt.hist(np.sum(np.where(a>threshold, a, 0), axis=1), bins=20, density=True)
#plt.hist(np.sum(np.where(a>threshold, a, 0), axis=0), bins=20, density=True)
#plt.hist(np.sum(np.where(a>threshold, a, 0), axis=0)+np.sum(np.where(a>threshold, a, 0), axis=1), bins=20, density=True)
plt.gca().set_xlabel('Strength', fontsize=20)
plt.xticks(range(0, len(matrices[-1])+1, 10))
plt.gca().set_ylabel('Population', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/distribution_strength.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Distribution of the strength', fontsize=25)
    plt.show()
    
    
#Distribution of the natural frequencies
plt.hist(naturalFrequencies, bins=30, density=True)
plt.gca().set_xlabel('$\omega$', fontsize=20)
plt.gca().set_ylabel('Population', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/natural_frequencies.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Distribution of the natural frequencies', fontsize=25)
    plt.show() 
       

    
#Distribution of the phase after updating
plt.hist(phases[-1,:], bins=20, density=True)
plt.gca().set_xlabel('$\\theta$', fontsize=20)
plt.xticks([0, pi, 2*pi], ['0', '$\pi$', '$2\pi$'])
plt.gca().set_ylabel('Population', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/distribution_phases.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Distribution of the phase', fontsize=25)
    plt.show() 
    

#The phase in function of the number of link
plt.plot(adjacency.sum(axis=1), phases[-1,:], 'o', color='red') 
plt.gca().set_ylim([0, 2*pi])
plt.gca().set_ylabel('Phase', fontsize=20)
plt.gca().set_xlabel('Degree', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/nodeDegree_Phases.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Relationship between the node degree and the phases', fontsize=25)
    plt.show() 
    
    
#Evolution of the autocorrelations
plt.plot(tpoints, autocorrelations) 
plt.gca().set_ylim([0, 1])
plt.gca().set_xlim([dt*len(phases)//5, dt*len(phases)])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
if save : 
    plt.savefig('results/autocorrelations.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Autocorrelations of phases', fontsize=25)
    plt.show() 


#Phases evolution
plt.matshow(np.transpose(phases), cmap=plt.cm.viridis, vmin=0, vmax=2*pi, extent=[0, int(len(phases)*dt), 0, len(phases[0])], origin='lower', aspect='auto')
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
plt.gca().set_ylabel('Neurons', fontsize=20)
plt.gca().xaxis.tick_bottom()
cbar = plt.colorbar(ticks=[0, pi, 2*pi])
cbar.ax.set_ylabel('$\\theta$', rotation=180)
cbar.ax.set_yticklabels(['0', '$\pi$', '$2\pi$'])  # vertically oriented colorbar
cbar.ax.tick_params(labelsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.gca().set_xlim(19000*dt, 20000*dt)
    plt.savefig('results/Phase_patterns_zoom_T0.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(60000*dt, 64000*dt)
    plt.savefig('results/Phase_patterns_zoom_T1.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(249000*dt, 250000*dt)
    plt.savefig('results/Phase_patterns_zoom_T2.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.savefig('results/Phase_patterns.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Phase patterns', fontsize=25)
    plt.show() 


#Phases evolution sorted
plt.matshow(np.transpose(phases)[order], cmap=plt.cm.viridis, vmin=0, vmax=2*pi, extent=[0, int(len(phases)*dt), 0, len(phases[0])], origin='lower', aspect='auto')
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
plt.gca().set_ylabel('Neurons', fontsize=20)
plt.gca().xaxis.tick_bottom()
cbar = plt.colorbar(ticks=[0, pi, 2*pi])
cbar.ax.set_ylabel('$\\theta$', rotation=180)
cbar.ax.set_yticklabels(['0', '$\pi$', '$2\pi$'])  # vertically oriented colorbar
cbar.ax.tick_params(labelsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.gca().set_xlim(299000*dt, 300000*dt)
    plt.savefig('results/Phase_patterns_zoom_T3.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.savefig('results/Phase_patterns_sorted.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Phase patterns (sorted)', fontsize=25)
    plt.show() 


 
#Limit cyle represented in mean phase space
plt.scatter(np.mean(phases, axis=1), np.roll(np.mean(phases, axis=1), -1), s=2, color='black')
plt.gca().set_ylabel('m(t+1)', fontsize=20)
plt.gca().set_xlabel('m(t)', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
   plt.savefig('results/limit_cyle.png', dpi=300, bbox_inches='tight')
   plt.close()
else :
   plt.gca().set_title('Limit cyle represented in mean phase space', fontsize=25)
   plt.show()
    

    
""" 
#Evolution of the average phases/potential through the time
plt.plot(tpoints, np.mean(phases, axis=1))
plt.gca().set_ylabel('Potential/Phase', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
if save : 
   plt.savefig('results/potential_system_evolution.png', dpi=300, bbox_inches='tight')
   plt.close()
else :
   plt.gca().set_title('Average potential of the network', fontsize=25)
   plt.show()
"""



# Mutual Information and entropy
phases_initial = phases[0,:]
phases_final = phases[-1,:]    

# Discretization
bins = 10
phases_initial_binned = np.digitize(phases_initial, np.linspace(0, 2 * np.pi, bins))
phases_final_binned = np.digitize(phases_final, np.linspace(0, 2 * np.pi, bins))

mi = mutual_info_score(phases_initial_binned, phases_final_binned)
print("Mutual Information:", mi)

prob, _ = np.histogram(phases_final_binned, bins=bins, density=True)
H = entropy(prob, base=2)  # base=2 pour obtenir des bits
print("Entropy of final phase pattern:", H)



"""
#Trajectory
delta_phi = phases[:,0] - phases[:,1]    

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

ax.plot(matrices[:,0,1], matrices[:,1,0], delta_phi, color='k', lw=2)

ax.set_xlabel('k12')
ax.set_ylabel('k21')
ax.set_zlabel('Δφ')

plt.show()
"""



"""
# Lyapunov_exponents
perturbation_size = 1e-6 
time_values = np.linspace(0, 1.0, len(phases))
num_lyapunov = len(phases[0])  
perturbations = np.eye(len(phases[0])) * perturbation_size 

lyapunov_exponents = np.zeros(num_lyapunov)

for t in range(1, len(phases)):
    delta_phases = phases[t] - phases[t - 1]
    
    coupling_matrix = np.outer(delta_phases, delta_phases) 
    
    new_perturbations = perturbations + coupling_matrix @ perturbations
    
    Q, R = np.linalg.qr(new_perturbations)
    perturbations = Q 
    
    lyapunov_exponents += np.log(np.abs(np.diag(R)))

lyapunov_exponents /= time_values[-1]

print("Lyapunov exponents:", lyapunov_exponents)
"""

