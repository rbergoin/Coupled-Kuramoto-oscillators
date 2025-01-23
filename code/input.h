#include "graph.h"



/********************** Additional processing ************************/


/*
* Function to swap elements of array
*
* @param	xp		array 1
* @param	yp		array 2
*/
void swap(double* xp, double* yp) 
{ 
    double temp = *xp; 
    *xp = *yp; 
    *yp = temp; 
} 


/*
* Function to save a spike of neuron
*
* @param	i			index of the neuron
* @param	t_spike		time of the spike in second
*/
void saveSpike(int i, long double t_spike)
{
	/**** Save spikes  through the time ****/
	FILE *fptr = fopen("spikes.txt","a");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	fprintf(fptr,"%d %3.5Lf\n", i, t_spike);
	
	fclose(fptr);
}
 
 
/*
* Function to perform Selection Sort 
*
* @param	arr		array to sort
* @param	n		number of element
*/
void selectionSort(double arr[], int n) 
{ 
	int i, j, min_idx; 
  
    //One by one move boundary of unsorted subarray 
    for (i = 0; i < n - 1; i++) 
	{ 
		//Find the minimum element in unsorted array 
		min_idx = i; 
		for (j = i + 1; j < n; j++)
		{
			if (arr[j] < arr[min_idx])
			{ 
				min_idx = j;
			} 
  		} 
        //Swap the found minimum element 
        //with the first element 
        swap(&arr[min_idx], &arr[i]); 
    } 
}




/********************** Filled vector ************************/


/*
* Change the input vector by null values
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
*/
void addNullInputs(double *inputs, int n)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		inputs[i] = 0.0;
	}
}


/*
* Change the input vector by binary values of one
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
*/
void addBinaryInputs(double *inputs, int n)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		inputs[i] = 1.0;
	}
}


/*
* Change the input vector by binary values of one
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
*/
void addBinaryRandom2(double *inputs, int n)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		inputs[i] = get_random(3.0, 4.0);
	}
}


/*
* Change the input vector by 2 values of 1 and 2
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
*/
void add2Inputs(double *inputs, int n, bool *inhibitory)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		/*
		if((i%2==0) && !inhibitory[i])
		{
			inputs[i] = 1.0;
		}
		else if(!inhibitory[i])
		{
			inputs[i] = 3.0;
		}
		*/
		
		if(i<(n/2) && !inhibitory[i])
		{
			inputs[i] = 1.0;
		}
		else if(!inhibitory[i])
		{
			inputs[i] = 3.0;
		}
	}
}


/*
* Change the input vector by values obtained with a Gaussian distribution, mean = 2, standard deviation = 0.2
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
*/
void addGaussianInputs(double *inputs, int n)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		inputs[i] = get_random_normal(2.0, 0.2, 1.0, 3.0);
	}
	
	selectionSort(inputs, n); 
}


/*
* Change the input vector by values uniformly distributed between 0,75 and 1,25
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
*/
void addUniformInputs(double *inputs, int n)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		inputs[i] = 1.0 + (1.0/(2.0*n))*(i-(n/2.0));
	}
}


/*
* Change the input vector by linear values from min to max located from start to end
*
* @param	inputs			the input vector
* @param	min				min value
* @param	max				max value
* @param    start			index of start of input
* @param    end				index(excluded) of end of input
*/
void addLinearInputsBounded(double *inputs, float min, float max, int start, int end)
{
	int i;
	double j = 0.0;
	
	//For all inputs
	for(i=start; i < end; i++)
	{
		inputs[i] = (j/(end-start-1.0))*(max-min)+min;
		j++;
	}
}




/********************** Binary values ************************/


/*
* Change the input vector by binary values on a specific cluster, stimulation brain 50-200hz
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	index			the index of the center position the cluster
* @param	nbInput			the number of input to activatee
* @param	tau_m			the membrane time constant of the target neurons
* @param	inputFrequency	the frequency of the external input applied
*/
void addBinaryLocalized(double *inputs, int n, int index, int nbInput, float tau_m, float inputFrequency)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		//If we are on the cluster
		if((i>=(index-nbInput/2)) && (i<(index+nbInput/2)))
		{
			inputs[i] = pow(inputFrequency*M_PI*tau_m, 2.0);	//inputFrequency Hz;
			//inputs[i] = get_random(pow(50.0*M_PI*tau_m, 2.0), pow(inputFrequency*M_PI*tau_m, 2.0));
			//inputs[i] = sin(t*tau_m)*sin(t*tau_m)*pow(50.0*M_PI*tau_m, 2.0);	//200 Hz;
			//freq = rand()%201;
			//inputs[i] = pow(freq*M_PI*tau_m, 2.0);	
		}
		else
		{
			inputs[i] = 0.0;
		}
	}
}


/*
* Change the input vector by binary values on a specific cluster, stimulation brain 50-200hz
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	index			the index of the center position the cluster
* @param	nbInput			the number of input to activate
* @param	inhibitory		if neuron is inhibitory or not	
*/
void addBinaryLocalizedNoReset(double *inputs, int n, int index, int nbInput, bool *inhibitory)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		//If we are on the cluster
		if((i>=(index-nbInput/2)) && (i<(index+nbInput/2)) && !inhibitory[i])
		{
			//inputs[i] = (50.0*2.0*M_PI)/100.0;
			inputs[i] = pow(50.0*(M_PI/100.0), 2.0);
			//inputs[i] = sqrt(3.0)*2.0;
			//inputs[i] = 3.0;
		}
	}
}


/*
* Change the input vector by null values on a specific cluster
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	index			the index of the center position the cluster
* @param	nbInput			the number of input to activate
* @param	inhibitory		if neuron is inhibitory or not	
*/
void addResetLocalized(double *inputs, int n, int index, int nbInput, bool *inhibitory)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		//If we are on the cluster
		if((i>=(index-nbInput/2)) && (i<(index+nbInput/2)) && !inhibitory[i])
		{
			inputs[i] = 0;
		}
	}
}





/*
* Change the input vector by binary values on a specific cluster with random activation
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	index			the index of the center position the cluster
* @param	nbInput			the number of input to activate
* @param	tau_m			the membrane time constant of the target neurons
* @param	inputFrequency	the frequency of the external input applied
*/
void addBinaryLocalizedRandom(double *inputs, int n, int index, int nbInput, float tau_m, float inputFrequency)
{
	int i;
	//int ratio = (rand()%2) +1;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		//If we are on the cluster
		if((i>=(index-nbInput/2)) && (i<(index+nbInput/2)))
		{
			if((rand()%2)!=0)
			//if((rand()%ratio)==0)
			{
				//inputs[i] = (50.0*2.0*M_PI)/100.0;
				//inputs[i] = pow(50.0*(M_PI/100.0), 2.0);
				//inputs[i] = sqrt(3.0)*2.0;
				//inputs[i] = 3.0;
				inputs[i] = pow(inputFrequency*M_PI*tau_m, 2.0);	//200 Hz;
			}
			else
			{
				inputs[i] = 0.0;
			}
		}
		else
		{
			inputs[i] = 0.0;
		}
	}
}





/*
* Change the input vector by binary values on a specific cluster and on hubs 
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	nbClusters		the number of clusters in the vector
* @param	index			the index of the cluster selected
* @param	nbInput			the number of input to activate
* @param	hubs			the number of hubs in the vector
*/
void addBinaryLocalizedHubs(double *inputs, int n, int index, int nbInput, int hubs)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		//If we are on the cluster
		if((i>=(index-nbInput/2)) && (i<(index+nbInput/2)))
		{
			inputs[i] = 1.0;
		}
		//If we are on the hubs
		else if(i>=(n-hubs))
		{
			inputs[i] = 1.0;
		}
		else
		{
			inputs[i] = 0.0;
		}
	}
}


/*
* Change the input vector by binary values on random positions according to an uniform distribution
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	nbInput			the number of input to activate randomly
*/
void addBinaryRandom(double *inputs, int n, int nbInput)
{
	int i, index;
	
	//Clear all inputs 
	for(i=0; i < n; i++)
	{
		inputs[i] = 0.0;
	}
	
	//For nbInput activated
	for(i=0; i < nbInput; i++)
	{
		do
		{
			index = rand()%n;
		} while(inputs[index] == 1.0) ;
		inputs[index] = 1.0;
	}
}


/*
* Change the input vector by binary values on random positions according to a gaussian distribution
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	nbInput			the number of input to activate randomly
*/
void addBinaryRandomGaussian(double *inputs, int n, int nbInput)
{
	int i, index;
	
	//Clear all inputs 
	for(i=0; i < n; i++)
	{
		inputs[i] = 0.0;
	}
	
	//For nbInput activated
	for(i=0; i < nbInput; i++)
	{
		do
		{
			index = (int)get_random_normal(n/2.0, n/4.0, 0, (n-1));
		} while(inputs[index] == 1.0) ;
		inputs[index] = 1.0;
	}
}




/********************** Analog values ************************/


/*
* Change the input vector by gaussian values on a specific cluster
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	mu 				the mean / position of the gaussian
* @param	sigma 			the stantard deviation of the gaussian
* @param	A				the amplitude of the gaussian
* @param	nbInput			the number of input to activate
*/
void addPartialGaussianLocalized(double *inputs, int n, float mu, float sigma, float A, int nbInput)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		//If we are on the cluster
		if((i>=(mu-nbInput/2)) && (i<=(mu+nbInput/2)))
		{
			inputs[i] = A * (1.0/(sigma*sqrt(2.0*M_PI))) * exp(-pow(i-mu, 2.0)/(2.0*sigma*sigma));
		}
		else
		{
			inputs[i] = 0.0;
		}
	}
}


/*
* Change the input vector by gaussian values on a specific cluster
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	mu 				the mean / position of the gaussian
* @param	sigma 			the stantard deviation of the gaussian
* @param	A				the amplitude of the gaussian
*/
void addGaussianLocalized(double *inputs, int n, float mu, float sigma, float A)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		inputs[i] = A * (1.0/(sigma*sqrt(2.0*M_PI))) * exp(-pow(i-mu, 2.0)/(2.0*sigma*sigma));	
	}
}


/*
* Change the input vector by values distributed with a gaussian in a localized area
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	index			the index of the center position the cluster
* @param	nbInput			the number of input to activate
*/
void addGaussianDistributionLocalized(double *inputs, int n, int index, int nbInput)
{
	int i, j = 0;
	double newInputs[nbInput];
	
	//For all new inputs
	for(i=0; i < nbInput; i++)
	{
		newInputs[i] = get_random_normal(2.0, 0.2, 1.0, 3.0);
	}
	selectionSort(newInputs, nbInput); 
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		//If we are on the cluster
		if((i>=(index-nbInput/2)) && (i<=(index+nbInput/2)))
		{
			inputs[i] = newInputs[j];
			j++;
		}
		else
		{
			inputs[i] = 0.0;
		}
	}
}


/*
* Change the input vector by values distributed uniformly in a localized area
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	index			the index of the center position the cluster
* @param	nbInput			the number of input to activate
*/
void addUniformDistributionLocalized(double *inputs, int n, int index, int nbInput)
{
	int i, j = 0;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		//If we are on the cluster
		if((i>=(index-nbInput/2)) && (i<=(index+nbInput/2)))
		{
			inputs[i] =  1.0 + (1.0/(2.0*nbInput))*(j-(nbInput/2.0));
			j++;
		}
		else
		{
			inputs[i] = 0.0;
		}
	}
}


