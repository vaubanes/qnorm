/******************************************************************************
* FILE: Qnorm.c
* DESCRIPTION:
*   Improved sequential prototype for Qnorm 
*   Qnormalisation Method: function that implements the ben Bolstad Method
*   quantile Normalization of High density Oliglonucleotide Array Data
*
* AUTHOR: O.Trelles (23 Feb.09)
* 23.Feb.09  : Using command line argums
*              Qnorm fMatrix.IN nRow nExp Normalised.fname  mode
*
*              fList  : contains a list of nExp filenames (gene-expression data files) 
*                          line format: fileName[TAB]nGenes[TAB]FileType[NEWLINE]
*              nRows  : number of genes in each file 
*              Normalised.fname where the normalised values will be stored 
*                       (as a text-tabulated matrix)
*              mode : m: keep the Index matrix in memory ( d: in disk)
*                       THIS VERSION REQUIRES ALL FILES WITH THE SAME NUMBER OF GENES
*
*
*  Command line Parameters
*         sintaxis:    Qnorm [-Option=value]... 
*
*  Option  Description              Default value     alternative Values
*  ------  -----------------------  -------------     ------------------
*  -p      NUmber of Nodes             4               not used
*  -i      File name (list of files)   qInput.txt      valid existing pathname
*  -o      Output binary matrix        qOut.bin        binary by columns file
*  -e      Number of experiments       2               positive integer
*  -g      NUmber of genes             15              positive integer
*  -t      Traspose the fileOut        Not             -T (yes)
*  -M      Index Matrix in mem         D (in disk)     -M (in memory)
*  -V      Verbose mode                Not             -V  
* ---------------------------------------------------------------------------
*
*	@returns >0 if everything was fine <0 if there was an error
* 22.04.09: syncho the pQnorm6 version with this sequential version
* LAST REVISED: 27/02/09
******************************************************************************/

#include "Qfunc.c"
#include "quicksort/quicksort.c"

int main(int ac, char **av){

        struct Files *fList=NULL;
        struct params *p=NULL;

        p = CommandLine(ac,av);

	if ((fList=LoadListOfFiles(p))==NULL) 
           terror("Loading list of files");

        QNormMain(p,fList);	

	
	return 1;
}


void QNormMain(struct params *p, struct Files* fList){

   double *dataIn, *dataOut;
   int **mIndex;
   int *dIndex;
   struct Average *AvG; // global Average by row
   int i,j,k;
   FILE *fI, *fOut;
   int nG=p->nG;
   int nE=p->nE;

   // Memory===========================================
   // Index array 
   if (p->MemIndex) { // in memory - full
     if ((mIndex=(int **)calloc(nG,sizeof(int*)))==NULL) 
        terror("memory for index1");
     for (i=0; i<nG;i++)
        if((mIndex[i]=(int *)calloc(nE,sizeof(int)))==NULL) 
           terror("memory for index2 full matrix");
   } else 
     if ((fI=fopen("~tmp","wb"))==NULL) 
        terror("opening tmp-index file");

   // This will always be necessary to decuple the function
   if((dIndex=(int *)calloc(nG,sizeof(int)))==NULL) 
       terror("memory for index2");


   if ((dataIn=(double *)calloc(nG,sizeof(double)))==NULL) 
       terror("memory for dataIn array");
   if ((AvG   =(struct Average *)calloc(nG,sizeof(struct Average)))==NULL) 
       terror("memory for Average array");

   for (j=0; j< nG;j++) { // init Accumulation array 
      AvG[j].Av=0;        // =HUGE_VAL; ???
      AvG[j].num=0;
   }


   // QNORM ===============================================================
   if (p->Verbose) fprintf(stderr,"kickoff..%d files.\n",nE);
   for (i=0; i< nE; i++) { // Qnorm for each datafile: STEP 1
        LoadFile(fList, i, dataIn);

        if (p->Verbose) fprintf(stderr,"%4d file: %s\n",i, fList[i].fname); 

        Qnorm1(dataIn, dIndex, fList[i].nG); // dataIn returns ordered and Index contains the origial position

        AccumulateRow(AvG, dataIn , nG);

        // now decide how to proceed with indexes
        if (p->MemIndex) { // in memory - full
          for (j=0;j<nG;j++) 
            mIndex[j][i]= dIndex[j];
        } else {         // in disk
          fseek(fI, nG*i*sizeof(int), SEEK_SET);
          fwrite(dIndex, sizeof(int), nG, fI);
        }


   }

   // Row average  ----------------------------------------------

   for (i=0;i<nG;i++) 
       AvG[i].Av /=AvG[i].num;




   // Finally produce the ORDERED output file [STEP 2]-------------------------

   if (!p->MemIndex) { 
      fclose(fI);
      if ((fI=fopen("~tmp","rb"))==NULL) terror("opening tmp-index for reading file");
   }

   // oputput file (by cols)
   if ((fOut=fopen(p->fOutName,"wb"))==NULL) terror("opening OUTPUT file");

   if ((dataOut=(double *)calloc(nG,sizeof(double)))==NULL) terror("memory for dataOut array");

   for (i=0;i<nE;i++) {
        if (p->MemIndex) { // in memory - full
          for (j=0;j<nG;j++) 
            dIndex[j]=mIndex[j][i];
        } else {
          fseek(fI, nG*i*sizeof(int), SEEK_SET);
          fread(dIndex, sizeof(int), nG, fI);
        }

        // complete the output vector
        for (j=0;j<nG;j++) 
          dataOut[dIndex[j]]=AvG[j].Av; // OJO

        fseek(fOut, (long)nG*i*sizeof(double), SEEK_SET);
        fwrite(dataOut, sizeof(double), nG, fOut);
           
   }
   fclose(fOut);
   if (!p->MemIndex) fclose(fI);

   return;

}


// input returns ordered and Index contains the origial position

int Qnorm1(double *input, int *dIndex, int nG){
	int i,j,k,n;

	for (j=0; j<nG;j++) dIndex[j]=j; // init the indexes array

// Quicksort (recursive) QsortC(input,0,nG-1,dIndex); 

        quicksort(input,dIndex,nG); // iterative quicksort (stack problems using recursivity)

        return 1;
}

void AccumulateRow(struct Average *AvG, double *input , int nG){
        int i;
 	
        for (i=0;i<nG;i++) {
		AvG[i].Av+=input[i];
		AvG[i].num++;
	}


	return;

}

