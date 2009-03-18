/******************************************************************************
* FILE: pQnorm1.c
* DESCRIPTION:
*   First parallel prototype for Qnorm using "Sections Work-sharing"
* AUTHOR: O.Trelles (23 Feb.09)
* 23.Feb.09  : Using command line argums
               pQnorm1 fMatrix.IN nRow nCol






* LAST REVISED: 23/02/09
******************************************************************************/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int Qnormalization(double **input, double **output, int nrows, int ncolumns);
void QsortC(double *array,int l,int r,int *index);
int partition( double* a, int l, int r, int *indexes);
void terror(char *);

void terror(char *s) {
  printf("\nERR**** %s ***\n");
  exit (-1);
}

/**
	Qnormalisation Method: function that implements the ben Bolstad Method
	quantile Normalization of High density Oliglonucleotide Array Data
	@input : data matrix (each column is an experiment, each row a gene or data)
	@output: data matrix with the results
	@nrows: number of genes or rows in the matrix
	@ncolumns: number of experiments in the matrix
	@returns >0 if everything was fine <0 if there was an error
**/
int Qnormalization(double **input, double **output, int nrows, int ncolumns)
{
	int i,j,k,n;
	double average;
	double** ordered;
	int** indexes;

	//initialize auxiliar array to not modify the input array
	ordered=(double **)malloc(nrows*sizeof(double*));
	for (i=0; i<nrows;i++){
		ordered[i]=(double *)malloc(ncolumns*(sizeof(double)));
	}
	//auxiliar array to store the indexes
	indexes=(int **)malloc(nrows*sizeof(int*));
	for (i=0; i<nrows;i++){
		indexes[i]=(int *)malloc(ncolumns*(sizeof(int)));
	}

	//copy the values 
	for (i=0; i<nrows;i++){	
		for (j=0; j<ncolumns;j++){
			ordered[i][j]=input[i][j];	
			//UNIFY NAN CONSTANT 
			if ((!__finite(input[i][j]))||(__isnan(input[i][j])))
				ordered[i][j]=HUGE_VAL;
		}
		
	}

	//initialize the indexes array
	for (i=0; i<nrows;i++){
		for (j=0; j<ncolumns;j++){
			indexes[i][j]=j;
		}
	}
	
	
	//order each data column
	for (i=0; i<nrows; i++){
		//Quicksort 
		QsortC(ordered[i],0,ncolumns-1,indexes[i]);
	//	bubbleSort(ordered[i], ncolumns, indexes[i]);
	}

	
	//APPLY QNORMALISATION
	//Take the means across rows of Xsort and assign this mean to each element in the row to get X0sort
    //. Get Xnormalized by rearranging each column of X0 sort to have the same ordering as original X	
	//calculate the average of each j elem
	for (j=0; j<ncolumns; j++){
		average=HUGE_VAL;
		n=0;
		for (i=0; i<nrows;i++){
			if ((__finite(ordered[i][j]))&&(!__isnan(ordered[i][j]))){
				if ((!__finite(average))||(__isnan(average))){
					average=0;
				}
				average+=ordered[i][j];
				n++;
			}
		}
		if ((__finite(average))&&(!__isnan(average))){
			average/=n;
		}
		for (i=0; i<nrows;i++){
			if ((__finite(ordered[i][j]))&&(!__isnan(ordered[i][j])))
				ordered[i][j]=average;
		}
	}
	
	//reorder each data column
	for (i=0; i<nrows;i++){
		for (j=0; j<ncolumns;j++){
			k=indexes[i][j];
			output[i][k]=ordered[i][j];
		}
	}

	return 1;
}

void QsortC(double *array,int l,int r,int *index)
{
	
   int j;
   if( l < r ) 
   {
   	// divide and conquer
       j = partition( array, l, r,index);
	 //  j=(l+r)/2;
       QsortC( array, l, j-1,index);
       QsortC( array, j+1, r,index);
   }
	
}





int partition( double* a, int l, int r, int *indexes) {
   int i=l;
   int j=r+1;
   int k;
   double t,pivot;
   pivot = a[l];
   //i = l; 
   //j = r+1;
		
   while( 1)
   {
	   do{
			++i; 
	   }while( a[i] <= pivot && i <= r );
	   do{
		   --j; 
	   }while( a[j] > pivot );
   		if( i >= j ) break;
   		t = a[i]; a[i] = a[j]; a[j] = t;
		k=indexes[i];indexes[i]=indexes[j];indexes[j]=k;
   }
   t = a[l]; a[l] = a[j]; a[j] = t;
   k=indexes[l];indexes[l]=indexes[j];indexes[j]=k;
   return j;
}


/** DUMMY TEST **/
 
int main(int ac, char **av){

	double** aux=(double **)malloc(sizeof(double *)*10);
	int i,j;

        if (ac!=4) terror("Usage pQnorm1 fMatrix.in  nRows  nCOls");

        /* it has non sense to load the full matrix in memory	for (i=0; i<10;i++){
		aux[i]=(double *)malloc(sizeof(double)*10);
	}
	
	for (i=0; i<10;i++){
		for (j=0; j<10; j++){
			aux[i][j]=(10-j)*30+i*20;
		}
	}
	
	Qnormalization(aux,aux,10,10); 
	for (i=0; i<10;i++){
		for (j=0; j<10; j++){
			printf("%f\t",aux[i][j]);
		}
		printf("\n");
	}
	
	return 1;
}
