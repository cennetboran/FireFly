

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <memory.h>

#define DUMP	1
#define MAX_FFA	1000
#define MAX_D	1000

using namespace std;

int D = 1000;			
int n = 20;			
int MaxGeneration;		
int NumEval;			
int Index[MAX_FFA];		

double ffa[MAX_FFA][MAX_D];	
double ffa_tmp[MAX_FFA][MAX_D]; 
double f[MAX_FFA];		
double I[MAX_FFA];		
double nbest[MAX_D];          
double lb[MAX_D];		
double ub[MAX_D];		

double alpha = 0.5;		
double betamin = 0.2;           
double gama = 1.0;		

double fbest;			

typedef double (*FunctionCallback)(double sol[MAX_D]);


double cost(double sol[MAX_D]);
double sphere(double sol[MAX_D]);


FunctionCallback function = &cost;


double alpha_new(double alpha, int NGen)
{
	double delta;			
	delta = 1.0-pow((pow(10.0, -4.0)/0.9), 1.0/(double) NGen);
	return (1-delta)*alpha;
}


void init_ffa()
{
	int i, j;
	double r;

	
	for (i=0;i<D;i++)
	{
		lb[i] = 0.0;
		ub[i] = 2.0;
	}

	for (i=0;i<n;i++)
	{
		for (j=0;j<D;j++)
		{
			r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			ffa[i][j]=r*(ub[j]-lb[j])+lb[j];
		}
		f[i] = 1.0;			
		I[i] = f[i];
	}
}


void sort_ffa()
{
	int i, j;

	
	for(i=0;i<n;i++)
		Index[i] = i;

	
	for(i=0;i<n-1;i++)
	{
		for(j=i+1;j<n;j++)
		{
			if(I[i] > I[j])
			{
				double z = I[i];	
				I[i] = I[j];
				I[j] = z;
				z = f[i];			
				f[i] = f[j];
				f[j] = z;
				int k = Index[i];	
				Index[i] = Index[j];
				Index[j] = k;
			}
		}
	}
}


void replace_ffa()
{
	int i, j;


	for(i=0;i<n;i++)
	{
		for(j=0;j<D;j++)
		{
			ffa_tmp[i][j] = ffa[i][j];
		}
	}

	
	for(i=0;i<n;i++)
	{
		for(j=0;j<D;j++)
		{
			ffa[i][j] = ffa_tmp[Index[i]][j];
		}
	}
}

void findlimits(int k)
{
	int i;

	for(i=0;i<D;i++)
	{
		if(ffa[k][i] < lb[i])
			ffa[k][i] = lb[i];
		if(ffa[k][i] > ub[i])
			ffa[k][i] = ub[i];
	}
}

void move_ffa()
{
	int i, j, k;
	double scale;
	double r, beta;

	for(i=0;i<n;i++)
	{
		scale = abs(ub[i]-lb[i]);
		for(j=0;j<n;j++)
		{
			r = 0.0;
			for(k=0;k<D;k++)
			{
				r += (ffa[i][k]-ffa[j][k])*(ffa[i][k]-ffa[j][k]);
			}
			r = sqrt(r);
			if(I[i] > I[j])	
			{
				double beta0 = 1.0;
				beta = (beta0-betamin)*exp(-gama*pow(r, 2.0))+betamin;
				for(k=0;k<D;k++)
				{
					r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
					double tmpf = alpha*(r-0.5)*scale;
					ffa[i][k] = ffa[i][k]*(1.0-beta)+ffa_tmp[j][k]*beta+tmpf;
				}
			}
		}
		findlimits(i);
	}
}

void dump_ffa(int gen)
{
	cout << "Dump at gen= " << gen << " best= " << fbest << endl;
}


void help()
{
	cout << "Syntax:" << endl;
	cout << "  Firefly [-h|-?] [-l] [-p] [-c] [-k] [-s] [-t]" << endl;
	cout << "    Parameters: -h|-? = command syntax" << endl;
	cout << "				 -n = number of fireflies" << endl;
	cout << "				 -d = problem dimension" << endl;
	cout << "				 -g = number of generations" << endl;
	cout << "				 -a = alpha parameter" << endl;
	cout << "				 -b = beta0 parameter" << endl;
	cout << "				 -c = gamma parameter" << endl;
}

int main(int argc, char* argv[])
{
        int i;
        int t = 1;	

        
         for(int i=1;i<argc;i++)
         {
            if((strncmp(argv[i], "-h", 2) == 0) || (strncmp(argv[i], "-?", 2) == 0))
            {
    		help();
    		return 0;
            }
            else if(strncmp(argv[i], "-n", 2) == 0)         
            {
    		n = atoi(&argv[i][2]);
            }
            else if(strncmp(argv[i], "-d", 2) == 0)		
            {
    		D = atoi(&argv[i][2]);
            }
            else if(strncmp(argv[i], "-g", 2) == 0)		
            {
    		MaxGeneration = atoi(&argv[i][2]);
            }
            else if(strncmp(argv[i], "-a", 2) == 0)		
            {
    		alpha = atof(&argv[i][2]);
            }
            else if(strncmp(argv[i], "-b", 2) == 0)		
            {
    		betamin = atof(&argv[i][2]);
            }
            else if(strncmp(argv[i], "-c", 2) == 0)		
            {
    		gama = atof(&argv[i][2]);
            }
            else
            {
    		cerr << "Fatal error: invalid parameter: " << argv[i] << endl;
    		return -1;
            }
        }

        
	srand(1);

	
	init_ffa();
#ifdef DUMP
	dump_ffa(t);
#endif

	while(t <= MaxGeneration)
	{
		
		alpha = alpha_new(alpha, MaxGeneration);

		
		for(i=0;i<n;i++)
		{
                        f[i] = function(ffa[i]);                        
			I[i] = f[i];					
		}

		
		sort_ffa();
		
		replace_ffa();

		
		for(i=0;i<D;i++)
			nbest[i] = ffa[0][i];
		fbest = I[0];

		
		move_ffa();
#ifdef DUMP
		dump_ffa(t);
#endif
		t++;
	}

	cout << "End of optimization: fbest = " << fbest << endl;

	return 0;
}


double cost(double* sol)
{
	double sum = 0.0;

	for(int i=0;i<D;i++)
		sum += (sol[i]-1)*(sol[i]-1);

	return sum;
}

double sphere(double* sol) {
	int j;
	double top = 0;
	for (j = 0; j < D; j++) {
		top = top + sol[j] * sol[j];
	}
	return top;
}
