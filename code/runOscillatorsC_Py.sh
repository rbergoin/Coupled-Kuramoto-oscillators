#!/bin/sh

#Compile and execute C program and run visualization with python program
#Developped by Raphael BERGOIN
#Run :  ./runOscillatorsC_Py.sh

echo "Compilation..."
gcc -W -Wall -o simulationOscillators simulationOscillators.c -lm




echo "Running..."
# Choose one of these executions
# Parameters : number of iterations | number of neurons | epsilon | alpha | beta | time step | adjacency policy | weight policy | frequency policy | phase policy | plasticity policy 

#./simulationOscillators 											#Default

#FIG.1
#./simulationOscillators 10000 2 0.05 0.1 -0.7 0.1 f r c r s			#Symmetric
#./simulationOscillators 10000 2 0.05 0.1 -0.3 0.1 f r c r s 			#Bistable
#./simulationOscillators 10000 2 0.05 0.1 0.3 0.1 f r c r s 			#Asymmetric
#./simulationOscillators 10000 2 0.05 0.1 0.7 0.1 f r c r s 			#Chaos

#FIG.2 and 3d
./simulationOscillators 10000 200 0.05 0.1 -0.9 0.1 f r c r s 			#Two-cluster state / Hebbian learning 
#./simulationOscillators 10000 200 0.05 0.15 -0.1 0.1 f r c r s 		#Coherent state	   / STDP learning
#./simulationOscillators 10000 200 0.05 0.125 0.6 0.1 f r c r s 		#Chaotic state	   / Anti-Hebbian learning / Asymmetric

#FIG.3 Chaos b
#./simulationOscillators 10000 2 0.05 0.1 0.55 0.1 f r c r s 				











echo "Plotting..."
python3 plotMatrix.py 0


