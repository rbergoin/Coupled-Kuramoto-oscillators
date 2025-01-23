#include "input.h"

/*
* Struct to represent a network of N oscillators 
*
*/
struct oscillators 
{
    int n;								//number of oscillators
    float epsilon;						//the dynamic of the oscillators
	float alpha;						//the phase delay
	float beta;							//the weights parameter
	float h;							//the stepsize of the RK method
	int **a;							//adjacency matrix
	double **k;							//coupling weights matrix
	double *omega;						//natural frequencies of the oscillators
	double *phi;						//phases of the oscillators
};


/*
* Copy values of a table in a new one
*
* @param	table	table to copy
* @param	n		size of the table
*/
double *copyTable(double *table, int n) 
{ 
	int i;
	double *copy = (double *) malloc(sizeof(double) * n);
	
	for(i=0; i<n; i++)
	{
		copy[i] = table[i];
	}
	
	return copy;
}

/*
* Copy values of a matrix in a new one
*
* @param	matrix	matrix to copy
* @param	n		size of the matrix
*/
double **copyMatrix(double **matrix, int n) 
{ 
	int i, j;
	double **copy = (double **) malloc(sizeof(double *) * n);

	for (i=0; i<n; i++) 
	{
		copy[i] = (double *) malloc(sizeof(double) * n);
		for (j=0; j<n; j++) 
		{
			copy[i][j] = matrix[i][j];
		}
	}
	
	return copy;
}


/*
* Copy values of a matrix in a new one
*
* @param	matrix	matrix to copy
* @param	n		size of the matrix
*/
void freeMatrix(double **matrix, int n)  
{ 
	int i;
	
	//Free matrix
	for (i=0; i<n; i++) 
	{
		free(matrix[i]);
	}
	free(matrix);
}


/*
* Initialize a network of oscillators
*
* @param	n			number of oscillators
* @param	epsilon		the dynamic of the oscillators
* @param	alpha		the phase delay
* @param	beta		the weights parameter
* @param	h 		 	the stepsize of the RK method
* @param	adjacencyPolicy		type of policy for the adjacency matrix creation
* @param	weightPolicy		type of policy for the coupling weights matrix creation
* @param	frequencyPolicy		type of policy for the natural frequencies array creation
* @param	phasePolicy			type of policy for the phases array creation
*/
struct oscillators initOscillators(int n, float epsilon, float alpha, float beta, float h, char* adjacencyPolicy, char weightPolicy, char frequencyPolicy, char phasePolicy)
{    	
	int i, j;
	struct oscillators oscillators;
	oscillators.n = n;
	oscillators.epsilon = epsilon;
	oscillators.alpha = alpha;
	oscillators.beta = beta;
	oscillators.h = h;
	
	//Allocate and initialize adjacency matrix
	oscillators.a = (int **) malloc(sizeof(int *) * n);
	
	for (i=0; i<n; i++) 
	{
		oscillators.a[i] = (int *) malloc(sizeof(int) * n);
		for (j=0; j<n; j++) 
		{
			oscillators.a[i][j] = 0;							//By default no connected
			
			if(strcmp(adjacencyPolicy, "f")==0)					//Fully connected, no recurrent links
			{
				if(i!=j)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "fr")==0)			//Fully connected with recurrent links (all)
			{
				oscillators.a[i][j] = 1;
			}
			else if(strcmp(adjacencyPolicy, "fe")==0)			//Fully connected except one neuron with no connection
			{
				if(i!=0)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "re")==0)			//Randomly uniform connected except one neuron with no connection
			{
				if(i!=0)
				{
					oscillators.a[i][j] = rand()%2;
				}
			}
			else if(strcmp(adjacencyPolicy, "ru")==0)			//Randomly uniform connected
			{
				oscillators.a[i][j] = rand()%2;
			}
			else if(strcmp(adjacencyPolicy, "r23")==0)			//Randomly 2/3 connected
			{
				if((rand()%3)%2 == 0)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "r13")==0)			//Randomly 1/3 connected
			{
				if((rand()%3)%2 == 1)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "r14")==0)			//Randomly 1/4 connected
			{
				if((rand()%4) == 0)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "r34")==0)			//Randomly 3/4 connected
			{
				if((rand()%4) != 0)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "r120")==0)			//Randomly 1/20 connected
			{
				if((rand()%20) == 0)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "o")==0)			//One to one connected (in index order)
			{
				if(i == j+1)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "t")==0)			//Two by two connected (in index order)
			{
				if((i == j+1) || (i == j-1))
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "tr")==0)			//Two by two and recurrently connected (in index order)
			{
				if((i == j+1) || (i == j-1) || (i == j))
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "c")==0)			//Two by two in circle connected (in index order)
			{
				if((i == j+1) || (i == j-1))
				{
					oscillators.a[i][j] = 1;
				}
				if(i == 0)
				{
					oscillators.a[i][n-1] = 1;
				}
				if(i == (n-1))
				{
					oscillators.a[i][0] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "oa")==0)			//One(0) to all(1-N) connected (connected in index order)
			{
				if(i!=0)
				{
					oscillators.a[i][0] = 1;
				}
				if((i == j+1) || (i == j-1))
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "oaf")==0)			//One(0) to all(1-N) connected with feedback (connected in index order)
			{
				if(i!=0)
				{
					oscillators.a[i][0] = 1;
				}
				if((i == j+1) || (i == j-1))
				{
					oscillators.a[i][j] = 1;
					if(j!=0)
					{
						oscillators.a[0][j] = 1;
					}
				}
			}
			else if(strcmp(adjacencyPolicy, "ao")==0)			//All(1-N) to one connected (0)
			{
				if((i==0) && (j!=0))
				{
					oscillators.a[i][j] = 1;
				}
				else if((j==0) && (i!=0))
				{
					oscillators.a[i][j] = 1;
				}
			}
		}
	}
	
	if(strcmp(adjacencyPolicy, "sw")==0)			
	{
		createSmallWorld(oscillators.a, n/4.0, 0.1, n);					//Small world connected, K=n/4, p=0.1
	}
	else if(strcmp(adjacencyPolicy, "sf")==0)
	{
		createScaleFree(oscillators.a, (0.05*n), 0.05*n, n);			//Scale free connected, K=m0, m0=0.05n
	}
	else if(strcmp(adjacencyPolicy, "swr")==0)			
	{
		createSmallWorldRandom(oscillators.a, n/4.0, 0.1, n);			//Small world random connected, K=n/4, p=0.1
	}
	else if(strcmp(adjacencyPolicy, "sfr")==0)
	{
		createScaleFreeRandom(oscillators.a, 0.05*n, n);	//Scale free random connected, m0=0.05n
	}
	else if(strcmp(adjacencyPolicy, "mod")==0)
	{
		createModular(oscillators.a, n/4.0, n);							//Modular connected, K=n/4
	}
	else if(strcmp(adjacencyPolicy, "ml")==0)
	{
		createMultiLayer(oscillators.a, 4, n);							//Multi layer connected, nb layer = 4
	}
	else if(strcmp(adjacencyPolicy, "mlf")==0)
	{
		createMultiLayerFeedback(oscillators.a, 4, n);					//Multi layer connected, nb layer = 4
	}
	else if(strcmp(adjacencyPolicy, "af")==0)
	{
		createAllToFull(oscillators.a, 10, n);							//All to full connected, nb input = 10
	}
	else if(strcmp(adjacencyPolicy, "res")==0)
	{
		createReservoir(oscillators.a, 6, 2, n);						//Reservoirconnected, nb input = 6, nb output = 2
	}
		
	
	//Allocate and initialize coupling weights array
	oscillators.k = (double **) malloc(sizeof(double *) * n);

	for (i=0; i<n; i++) 
	{
		oscillators.k[i] = (double *) malloc(sizeof(double) * n);
		for (j=0; j<n; j++) 
		{
			if(oscillators.a[i][j]==1)  //if the link between node i and j exists
			{
				switch(weightPolicy)
				{
				case 'r' :	//Random value uniform between -1 and 1
					oscillators.k[i][j] = get_random(-1, 1);
					break;
				case 'p' :	//Constant value of +1
					oscillators.k[i][j] = 1.0;
					break;
				case 'm' :	//Constant value of -1
					oscillators.k[i][j] = -1.0;
					break;
				case 'n' :	//Constant value null
					oscillators.k[i][j] = 0.0;
					break;
				case 't' :	//Two values -1 +1
					if((i%2)==0)
					{
						oscillators.k[i][j] = -1.0;
					}
					else
					{
						oscillators.k[i][j] = 1.0;
					}
					break;
				default : //Random value uniform between -1 and 1
					oscillators.k[i][j] = get_random(-1, 1);
					break;
				}
			}
			else
			{
				oscillators.k[i][j] = 0;
			}
		}
	}	
	
	
	//Allocate and initialize natural frequencies array
	oscillators.omega = (double *) malloc(sizeof(double) * n);
	for(i=0; i<n; i++)
	{
		switch(frequencyPolicy)
		{
			case 'c' :	//Constant value of 1.0
				oscillators.omega[i] = 1.0;
				break;
			case 'e' :	//constant value of 1 except one exception with frequency of 3
				if(i!=0)
				{
					oscillators.omega[i] = 1.0;
				}
				else
				{
					oscillators.omega[i] = 3.0;
				}
				break;
			case 'r' :	//Random value between 1 and 2
				oscillators.omega[i] = get_random(1.0, 2.0);
				break;
			case 'n' :	//Random value between 0.0 and 2.0 distributed with a normal/gaussian law
				oscillators.omega[i] = get_random_normal(2.0, 0.2, 1.0, 3.0);
				break;
			case 'l' :	//Linear value close to w0=1.0, info = 1/2n
				oscillators.omega[i] = 1.0 + (1.0/(2.0*n))*(i-(n/2.0));
				break;
			case 'p' :	//Pow values between w0=1.0 and 2.0, info = 1/(n*n)
				oscillators.omega[i] = 1.0 + (1.0/(n*n))*i*i;
				break;
			case '2' :	//Two frequencies (1 and 2)
				oscillators.omega[i] = rand()%2 +1;
				break;
			case '3' :	//Three frequencies (1, 2, 3)
				oscillators.omega[i] = rand()%3 +1;
				break;
			case '4' :	//Four frequencies (1, 2, 3, 4)
				oscillators.omega[i] = rand()%4 +1;
				break;
			default : //Constant value of 1
				oscillators.omega[i] = 1.0;
				break;
		}
	}
	
	
	//Allocate and initialize phases array
	oscillators.phi = (double *) malloc(sizeof(double) * n);
	for(i=0; i<n; i++)
	{
		switch(phasePolicy)
		{
			case 'r' :	//Random value between 0 and 2PI
				oscillators.phi[i] = get_random(0, 2*M_PI); 
				break;
			case 'c' :	//Constant value of 0
				oscillators.phi[i] = 0.0;
				break;
			case '2' :	//Two phases (0 and PI)
				oscillators.phi[i] = (rand()%2)*M_PI;
				break;
			case '3' :	//Three phases (0, PI/2, PI)
				oscillators.phi[i] = (rand()%3)*M_PI/2.0;
				break;
			case 'o' :	//n phases ordered
				oscillators.phi[i] = ((double)i/n)*2*M_PI;
				break;
			default :	//Random value between 0 and 2PI
				oscillators.phi[i] = get_random(0, 2*M_PI); 
				break;
		}
	}
	
	return oscillators;
}



/*
* Free the memory of the oscillators
*
* @param	oscillators			pointer on the current network
*/
void freeOscillators(struct oscillators *oscillators)
{    	
	int i;
	
	//Free coupling weights array
	for (i=0; i<oscillators->n; i++) 
	{
		free(oscillators->k[i]);
	}
	free(oscillators->k);
	
	//Free adjacency matrix
	for (i=0; i<oscillators->n; i++) 
	{
		free(oscillators->a[i]);
	}
	free(oscillators->a);
	
	//Free natural frequencies array
	free(oscillators->omega);
	
	//Free phases array
	free(oscillators->phi);
}


/*
* Coupling function
*
* @param	k			coupling weight
* @param	phi			the phase
* @param	alpha		the phase delay
*/
double gammaF(double k, double phi, float alpha)
{
	return -k * sin(phi + alpha);
}
    
    
/*
* Coupling weights/plasticity function
*
* @param	phi			the phase
* @param	beta		the weights parameter
* @param	plasticityPolicy	type of policy for the plasticity function
*/   
double lambdaF(double phi, float beta, char plasticityPolicy)
{
	switch(plasticityPolicy)
	{
	case 's' :	//sinus plasticity function
		return -sin(phi + beta);
		break;
	case '3' :	//cosinus plasticity function with 3 cluster
		return cos(2*phi + beta);
		break;
	case '4' :	//cosinus plasticity function with 4 cluster
		return cos(3*phi + beta);
		break;
	default :	//sinus plasticity function
		return -sin(phi + beta);
		break;
	}
}


/*
* Runge-kutta method 4th order for the phase (phi)
*
* @param	phi_i		the phase for the oscillator i
* @param	phi_j		the phase vector for the neighbor of i
* @param	h			the stepsize of the RK method
* @param	a			the adjacency vector for i and its neighbors
* @param	k			the couping vector for i and its neighbors
* @param	omega		the nature frequency for the oscillator i
* @param	alpha		the phase delay
* @param	n			number of neighbors
*/ 
double RK_phi(double phi_i, double *phi_j, float h, int *a, double *k, double omega, float alpha, int n, int degree)
{   
	int j;
	double k1, k2, k3, k4;
	
	if(degree==0)
	{
		return (phi_i + h*omega);
	}
	
	double sum = 0.0;
	//for each neighbors of oscillator i
	for(j=0; j < n; j++)
    {
		sum += a[j]*gammaF(k[j], (phi_i - phi_j[j]), alpha);
	}
	k1 = omega + 1.0/degree * sum;
	
	sum = 0.0;
	//for each neighbors of oscillator i
	for(j=0; j < n; j++)
    {
		sum += a[j]*gammaF(k[j], (phi_i+k1*h/2.0 - phi_j[j]), alpha);
	}
	k2 = omega + 1.0/degree * sum;    

	sum = 0.0;
	//for each neighbors of oscillator i
	for(j=0; j < n; j++)
    {
		sum += a[j]*gammaF(k[j], (phi_i+k2*h/2.0 - phi_j[j]), alpha);
	}
	k3 = omega + 1.0/degree * sum;           
    
	sum = 0.0;
	//for each neighbors of oscillator i
	for(j=0; j < n; j++)
    {
		sum += a[j]*gammaF(k[j], (phi_i+k3*h - phi_j[j]), alpha);
	}
	k4 = omega + 1.0/degree * sum;
  
    return (phi_i + h*((k1 + 2*k2 + 2*k3 + k4)/6.0));
}    


/*
* Update oscillators' phases of the network
*
* @param	oscillators			pointer on the current network
*/   
void update_phases(struct oscillators *oscillators)
{
    int i;
	
	double *phases = copyTable(oscillators->phi, oscillators->n);	//copy phases before updating
	
	//for each oscillators of the network
	for(i=0; i < oscillators->n; i++)
    {   
		oscillators->phi[i] = RK_phi(oscillators->phi[i], phases, oscillators->h, oscillators->a[i], oscillators->k[i], oscillators->omega[i], oscillators->alpha, oscillators->n, nodeDegree(i, oscillators->a, oscillators->n));
	   
	   	//Limiting values
		if(oscillators->phi[i] >= (2.0*M_PI))	//phases must be between 0 and 2PI (periodic)
		{
			oscillators->phi[i] = oscillators->phi[i] - (2.0*M_PI);
		}
		else if(oscillators->phi[i] < 0)
		{
			oscillators->phi[i] = 0.0;
		}
    }
	
	free(phases); 	//free memory 
}





/*
* Update coupling weights of the network
*
* @param	oscillators			pointer on the current network
* @param	plasticityPolicy	type of policy for the plasticity function
*/  
void update_weights(struct oscillators *oscillators, char plasticityPolicy)
{
    int i, j;
	//for each oscillators of the network
	for(i=0; i < oscillators->n; i++)
    {
		for(j=0; j < oscillators->n; j++)
	    {
			if(oscillators->a[i][j]==1)  //if the link between from node j to i exists
			{
				oscillators->k[i][j] += oscillators->h * oscillators->epsilon * lambdaF((oscillators->phi[i] - oscillators->phi[j]), oscillators->beta, plasticityPolicy);	
				
				//Limiting values
				if(oscillators->k[i][j] > 1.0)
				{
					oscillators->k[i][j] = 1.0;
				}
				if(oscillators->k[i][j] < -1.0)
				{
					oscillators->k[i][j] = -1.0;
				}
			}
		}
	}
}


/*
* Change the natural frequencies of the network
*
* @param	oscillators			pointer on the current network
* @param	frequencies			new natural frequencies
* @param	n					number of frequencies to change
*/   
void set_natural_frequencies(struct oscillators *oscillators, double *frequencies, int n)
{
    int i;
	//for each oscillators of the network
	for(i=0; i < n; i++)
    {
		oscillators->omega[i] = frequencies[i];
	}
}