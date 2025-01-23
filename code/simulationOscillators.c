#include "oscillators.h"


/*
* Calculate the rate of change of weights
*
* @param	k		coupling weights matrix at time 1
* @param	k2		coupling weights matrix at time 2
* @param	n		number of oscillators
* @param	d		number of links
*/
double ChangeRate(double **k, double ** k2, int n, int d)
{
    int i, j;
    double K = 0.0;

	//for each oscillators of the network
	for (i=0; i<n; i++) 
	{
		for (j=0; j<n; j++) 
		{
			K += fabs(k2[i][j]-k[i][j]);
		}
	}
    
    return K/d;
}


/*
* Save a matrix in a file
*
* @param	fptr	file
* @param	m		matrix
* @param	l		number of line
* @param	c		number of column
*/
void saveMatrix(FILE *fptr, double **m, int l, int c)
{
	int i, j;
	
	//for each element of the matrix
	for(i=0; i < l; i++)
    {
		for(j=0; j < c; j++)
	    {
			fprintf(fptr,"%lf ", m[i][j]);
		}
		fprintf(fptr, "\n");
	}
	fprintf(fptr, "\n");
}


int main(int argc, char *argv[])
{	
    srand(time(NULL)); // randomize seed
	
	int t, i, j, degree;
	FILE *fptr;
	FILE *fptr2;
	FILE *fptr3;
	double **k1;
	int nbIterations;
	int n;					//number of oscillators
	float epsilon;			//the dynamic of the oscillators
	float alpha;			//the phase delay
	float beta;				//the weights parameter
	float h;				//the stepsize of the RK method	
	char* adjacencyPolicy;	//Type of policy for the adjacency matrix creation
	char weightPolicy;		//Type of policy for the coupling weights matrix creation
	char frequencyPolicy;	//Type of policy for the natural frequencies array creation
	char phasePolicy;		//Type of policy for the phases array creation
	char plasticityPolicy;	//Type of policy for the plasticity function
	
	if(argc < 12  || argc > 12) 
	{
		nbIterations = 1000;
		n = 200;
		epsilon = 0.005;
		alpha = 0.0;
		beta = 0.0;
		h = 1.0;
		adjacencyPolicy = "f";
		weightPolicy = 'r';
		frequencyPolicy = 'c';
		phasePolicy = 'r';
		plasticityPolicy = 's';
	}
	else
	{
		nbIterations = atoi(argv[1]);
		n = atoi(argv[2]);
		epsilon = atof(argv[3]);
		alpha = atof(argv[4])*M_PI;
		beta = atof(argv[5])*M_PI;
		h = atof(argv[6]);
		adjacencyPolicy = argv[7];
		weightPolicy = argv[8][0];
		frequencyPolicy = argv[9][0];
		phasePolicy = argv[10][0];
		plasticityPolicy = argv[11][0];
	}
	
	
	struct oscillators oscillators = initOscillators(n, epsilon, alpha, beta, h, adjacencyPolicy, weightPolicy, frequencyPolicy, phasePolicy); 		//Create a network of n oscillators
	
	degree = graphDegree(oscillators.a, n);
	

	
	/**** Save phases vector through the time ****/
	fptr = fopen("phases.txt","w");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	//for each oscillators of the network
	for(i=0; i < oscillators.n; i++)
    {
		fprintf(fptr,"%lf ",oscillators.phi[i]);
	}
	fprintf(fptr, "\n");
	
	
	/**** Save change rate of weights through the time ****/
	fptr2 = fopen("changeRates.txt","w");
	
	if(fptr2 == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	
	/**** Save weights matrix ****/
	fptr3 = fopen("weights_matrices.txt" ,"w");
	
	if(fptr3 == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	saveMatrix(fptr3, oscillators.k, oscillators.n, oscillators.n);
	
	
	
	/**** Simulation for nbIterations ****/
	for(t=0; t < nbIterations; t++)
	{	
		k1 = copyMatrix(oscillators.k, oscillators.n);	//copy weights before updating
		
		update_phases(&oscillators);

		update_weights(&oscillators, plasticityPolicy);
	
		
		
		/***** Save data *****/
		
		for(i=0; i < oscillators.n; i++)
	    {
			fprintf(fptr,"%lf ", oscillators.phi[i]);
		}
		fprintf(fptr, "\n");
		
		fprintf(fptr2,"%lf ", ChangeRate(k1, oscillators.k, oscillators.n, degree)); //calculate and register change rate of weights
		freeMatrix(k1, oscillators.n);
		
		if(((t+1)%1000)==0)	//save weight matrix every 1000 iterations
		{
			saveMatrix(fptr3, oscillators.k, oscillators.n, oscillators.n);
		}
	}
	
	fclose(fptr);
	fclose(fptr2);
	fclose(fptr3);
	
	
	/**** Save adjacency matrix ****/
	fptr = fopen("adjacency.txt" ,"w");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	//for each element of the matrix
	for(i=0; i < oscillators.n; i++)
    {
		for(j=0; j < oscillators.n; j++)
	    {
			fprintf(fptr,"%d ", oscillators.a[i][j]);
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);
	
	
	/**** Save natural frequencies ****/
	fptr = fopen("natural_frequencies.txt","w");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	//for each oscillators of the network
	for(i=0; i < oscillators.n; i++)
    {
		fprintf(fptr,"%f ", oscillators.omega[i]);
	}
	
	fclose(fptr);
	
	
	
	/**** Free memory ****/
	freeOscillators(&oscillators);
    
	
    return 0;
}