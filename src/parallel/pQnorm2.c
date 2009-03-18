/******************************************************************************
* FILE: pQnorm.c
* DESCRIPTION:
*   parallel version
*   Improved sequential prototype for Qnorm 
*   Qnormalisation Method: function that implements the ben Bolstad Method
*   quantile Normalization of High density Oliglonucleotide Array Data
*
* AUTHOR: O.Trelles (23 Feb.09)
* 23.Feb.09  : Using command line argums
*              Qnorm [-o=Value]  see below Command line params
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
*  -p      Number of Processors        4               Only in parallel version
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
*
* LAST REVISED: 27/02/09
******************************************************************************/
#include <omp.h>
#include "Qfunc.c"

int main(int ac, char **av){

        struct Files *fList=NULL;
        struct params *p=NULL;
        p = CommandLine(ac,av);

	if ((fList=LoadListOfFiles(p))==NULL) 
           terror("Loading list of files");
/*
        QNormMain(p,fList);	

	
	return 1;
}


void QNormMain(struct params *p, struct Files* fList){
*/
   double *dataIn, *dataOut;
   int **mIndex;
   int *dIndex;
   struct Average **AvG; // global Average by row (now a matrx)
   int i,j,k;
   FILE *fI, *fOut;
   int nG=p->nG;
   int nE=p->nE;
   int nP, nP1, tid;
   int From, To, Range,pos; 
   double checkVal=0;

// the partial AvG array used by each thread to get the 
// accumulated value, will be shared using a matrix (nPxAvGsize)
// to facilitate final global average

#pragma omp parallel shared(nG, nE,mIndex,AvG, fList, p) private(i,j,k,tid,nP1,dataIn, dataOut, dIndex,fOut,fI, From, To, Range)
 { // Open General parallel section [0]
    tid = omp_get_thread_num();
    nP = omp_get_num_threads();
    nP1=nP+1;
    
    if (p->Verbose) fprintf(stderr,"Thread %d starting...\n",tid);
    fprintf(stderr,"Thread tid=%d nP=%d\n",tid,nP);fflush(stderr);
#pragma omp sections nowait
    { // open Sections [1]
    #pragma omp section
    { // open Subsection [2a]
    if (tid == 0) {
       fprintf(stderr,"Aqui solo tid0=%d\n",tid);fflush(stderr);
       if (p->Verbose) fprintf(stderr,"Number of threads = %d\n", nP);
       // Memory===========================================
       // Only in the master?????
       // Index array 
       if (p->MemIndex) { // in memory - full
           if ((mIndex=(int **)calloc(nG,sizeof(int*)))==NULL) 
              terror("memory for index1");
           for (i=0; i<nG;i++)
             if((mIndex[i]=(int *)calloc(Range,sizeof(int)))==NULL) 
               terror("memory for index2 full matrix");
       }
       AvG =(struct Average **)calloc(nP1,sizeof(struct Average*));
       if (AvG==NULL) terror("memory for Average array");
       for (i=0; i<nP1;i++)
         if((AvG[i]=(struct Average *)calloc(nG,sizeof(struct Average)))==NULL) 
               terror("memory for average 2 full matrix");
       if ((fI=fopen("~tmp","wb"))==NULL) terror("opening tmp-index file");
       fprintf(stderr,"fin solo tid0=%d\n",tid);fflush(stderr);
     } // end tid=0
     fprintf(stderr,"aqui todos tid=%d\n",tid);fflush(stderr);

    } // close Subsection [2a]

    #pragma omp section
    { // open Subsection [2b]
   // LOAD DISTRIBUTION------------------
   Range = nE / nP;
   From = tid * Range;
   To   = (tid+1)*Range;
   if (To > nE) To = nE;

fprintf(stderr, "P5.. tid=%d\n",tid);fflush(stderr);

   // This will always be necessary to decuple the function
   if((dIndex=(int *)calloc(nG,sizeof(int)))==NULL) terror("memory for index2");

   if ((dataIn=(double *)calloc(nG,sizeof(double)))==NULL) terror("memory for dataIn array");

fprintf(stderr, "P5.5 tid = %d\n",tid);fflush(stderr);
   // each P initialize
   for (j=0; j< nG;j++) { // init Accumulation array 
      AvG[tid][j].Av=0;        // =HUGE_VAL; ???
      AvG[tid][j].num=0;
      if (tid==0) {
        AvG[nP][j].Av=0;        // =HUGE_VAL; ???
        AvG[nP][j].num=0;
      }
   }

fprintf(stderr, "P6\n");fflush(stderr);

   // QNORM ===============================================================
   if (p->Verbose) fprintf(stderr,"[1st-Step]");

   for (i=From; i< To; i++) { // Qnorm for each datafile: STEP 1
        LoadFile(fList, i, dataIn);
        if (p->Verbose) { fprintf(stderr,"."); fflush(stderr);}
        fprintf(stderr,"P8 tid=%d procesa exp i=%d\n",tid,i); fflush(stderr);

#ifdef DEBUG
        DebugPrint("Load", dataIn, fList[i].nG); 
#endif
        Qnorm1(dataIn, dIndex, fList[i].nG); // dataIn returns ordered and Index contains the origial position

#ifdef DEBUG
        DebugPrint("Sorted", dataIn, nG);         
#endif
        
        AccumulateRow(AvG[tid], dataIn , nG);

        // now decide how to proceed with indexes
        if (p->MemIndex) { // in memory - full
          for (j=0;j<nG;j++) 
            mIndex[j][i]= dIndex[j];
        } else {         // in disk
          pos=fList[i].pos;
          fseek(fI, nG*pos*sizeof(int), SEEK_SET);
          fwrite(dIndex, sizeof(int), nG, fI);
        }


#ifdef DEBUG
        fprintf(stderr,"Index (col=%d)\n",i);
        for (j=0;j<nG;j++) fprintf (stderr,"%d ", dIndex[j]); fprintf(stderr,"\n");
#endif

   } // end "for i" Qnorm for each datafile: STEP 1

   } // close Subsection [2b]

    #pragma omp section
    { // open subSection [2c]
      // Row average     [use col=tid=0] ------------------------------

fprintf(stderr, "P7 tid=%d\n",tid);fflush(stderr);
   if (tid==0) { // total in nP 
     for (i=0;i<nG;i++) {
       for (j=0;j<nP;j++){
/*       if ((__finite(AvG[i].Av))&&(!__isnan(AvG[i].Av))) 
           AvG[i].Av /=AvG[i].num; 
           else AvG[i].Av =0; 
*/
         AvG[nP][i].Av +=AvG[j][i].Av;
         AvG[nP][i].num+=AvG[j][i].num;
       }
       AvG[nP][i].Av /=AvG[nP][i].num;
       checkVal +=AvG[nP][i].Av;
      }
      checkVal /=nG;
      if (p->Verbose) fprintf(stderr, "checkVal=%lf\n",checkVal);
      fprintf(stderr, "tid=%d checkVal=%lf\n",tid,checkVal); fflush(stderr);

      // Now copy the global averge into the local ones
      // *distribute* the global array
      for (i=0;i<nG;i++) {
        for (j=0;j<nP;j++) {
           AvG[j][i].Av =AvG[nP][i].Av;
           AvG[j][i].num=AvG[nP][i].num;
        }
      }

     } // End tid=0


#ifdef DEBUG
     fprintf(stderr, "Row Average------------\n");
     for (j=0;j<nG;j++) fprintf (stderr,"%f (%d) ", AvG[j].Av,AvG[j].num); 
     fprintf(stderr,"\n");
#endif

fprintf(stderr, "P9 END tid=%d\n",tid);fflush(stderr);
  
  } // close subSection [2c]

 
  #pragma omp section
  { // open subSection [2d]
   // Finally produce the ORDERED output file [STEP 2]-------------------------

   if (!p->MemIndex) { 
      fclose(fI);
      if ((fI=fopen("~tmp","rb"))==NULL) terror("opening tmp-index for reading file");
   }

   // oputput file (by cols)
   if ((fOut=fopen(p->fOutName,"r+b"))==NULL) terror("opening OUTPUT file");

   if ((dataOut=(double *)calloc(nG,sizeof(double)))==NULL) terror("memory for dataOut array");
   if (p->Verbose) { fprintf(stderr,"\n[2nd Step]"); fflush(stderr);}

   for (i=From;i<To;i++) {
        if (p->Verbose) {fprintf(stderr,"."); fflush(stderr);}

        if (p->MemIndex) { // in memory - full
          for (j=0;j<nG;j++) 
            dIndex[j]=mIndex[j][i];
        } else {
          pos = fList[i].pos;
          fseek(fI, nG*pos*sizeof(int), SEEK_SET);
          fread(dIndex, sizeof(int), nG, fI);
        }

#ifdef DEBUG
        fprintf(stderr,"recovered Index (col=%d)\n",i);
        for (j=0;j<nG;j++) fprintf (stderr,"%d ", dIndex[j]); fprintf(stderr,"\n");
#endif


        // complete the output vector
        for (j=0;j<nG;j++) 
          dataOut[dIndex[j]]=AvG[tid][j].Av; // OJO

#ifdef DEBUG
        fprintf(stderr,"[pos=%ld] ",(long)nG*i*sizeof(double));
        DebugPrint("Out to write", dataOut, nG); 
#endif

        pos = fList[i].pos;
        fseek(fOut, (long)nG*pos*sizeof(double), SEEK_SET);
        fwrite(dataOut, sizeof(double), nG, fOut);
           
   }
   fclose(fOut);
   if (!p->MemIndex) fclose(fI);

   if (p->Traspose){
    if (p->Verbose) 
       { fprintf(stderr,"\nTransposing Disk->Disk"); fflush(stderr); }
    TransposeBin2Txt(p);
   }

   if (p->Verbose) fprintf(stderr,"Thread %d done.\n",tid);   

   } // close subSection [2d]

    } // close set of Sections [1]

 } // close General parallel section [0]
    return 1;

} // end MAIN------------------


// input returns ordered and Index contains the origial position

int Qnorm1(double *input, int *dIndex, int nG){
	int i,j,k,n;

	for (j=0; j<nG;j++) dIndex[j]=j; // init the indexes array

/*
	for (j=0; j<nG;j++) // UNIFY NAN CONSTANT 
	   if ((!__finite(input[j]))||(__isnan(input[j]))) input[j]=HUGE_VAL;
*/

	QsortC(input,0,nG-1,dIndex); // Quicksort 

        return 1;
}

void AccumulateRow(struct Average *AvG, double *input , int nG){
        int i;
 	
        for (i=0;i<nG;i++) {

/*
           if ((__finite(input[i])&&(!__isnan(input[i])))){
		if ((!__finite(AvG[i].Av))||(__isnan(AvG[i].Av))){
			AvG[i].Av=0;
		}
*/
		AvG[i].Av+=input[i];
		AvG[i].num++;
/*           }
*/	}


	return;

}

