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

   pQnorm3  : since the memory allocation is performed only by thread 0
              the starting point of the parallel section has been moved
              to the init of the main process----
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
   double *dataOut;
   int **mIndex;
   int *dInde2;
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

    tid = omp_get_thread_num();
    nP = p->nP;  // omp_get_num_threads();
    nP1=nP+1;
    
    if (p->Verbose) fprintf(stderr,"Thread %d starting...\n",tid);
    if (p->Verbose) fprintf(stderr,"Number of threads = %d\n", nP);
    fprintf(stderr,"Thread tid=%d nP=%d\n",tid,nP);fflush(stderr);
    // Memory===========================================
    // Only in the master?????
    // Index array 
    if (p->MemIndex) { // in memory - full
         if ((mIndex=(int **)calloc(nG,sizeof(int*)))==NULL) 
              terror("memory for index1");
         for (i=0; i<nG;i++)
           if((mIndex[i]=(int *)calloc(nE,sizeof(int)))==NULL) 
             terror("memory for index2 full matrix");
    } else {
        // open tmp file??
       if ((fI=fopen("~tmp","wb"))==NULL) terror("opening tmp-index file");
    }
    if ((AvG =(struct Average **)calloc(nP1,sizeof(struct Average*)))==NULL)
       terror("memory for Average array");
    for (i=0; i<nP1;i++)
       if((AvG[i]=(struct Average *)calloc(nG,sizeof(struct Average)))==NULL) 
              terror("memory for average 2 full matrix");
    fprintf(stderr,"fin solo tid0=%d\n",tid);fflush(stderr);


#pragma omp parallel shared(nG, nE,mIndex,AvG, fList, p) private(i,j,k,tid,nP1, dataOut,fOut,fI, From, To, Range)
 { // Open General parallel section [0]

    tid = omp_get_thread_num();
    if (!tid) {
      nP = omp_get_num_threads();
      if (nP != p->nP) terror("something wrong in nP");
    }

fprintf(stderr, "P4. tid=%d hay %d thread\n",tid,nP);fflush(stderr);
#pragma omp sections
//#pragma omp sections nowait
  { // open Sections [1]
 #pragma omp section 
    { // open Subsection [2b]

   double *dataIn;
   int *dIndex;

printf( "P4.4tid=%d en parallel\n",tid);
   // LOAD DISTRIBUTION------------------
   Range = nE / nP;
   From = tid * Range;
   To   = (tid+1)*Range;
   if (To > nE) To = nE;


   // This will always be necessary to decuple the function
   if((dIndex=(int *)calloc(nG,sizeof(int)))==NULL) terror("memory for index2");

   if ((dataIn=(double *)calloc(nG,sizeof(double)))==NULL) terror("memory for dataIn array");

fprintf(stderr, "P5. tid=%d fromTo=%d/%d\n",tid,From,To);fflush(stderr);
   // each P initialize
   for (j=0; j< nG;j++) { // init Accumulation array 
      AvG[tid][j].Av=0;        // =HUGE_VAL; ???
      AvG[tid][j].num=0;
      if (tid==0) {
        AvG[nP][j].Av=0;        // =HUGE_VAL; ???
        AvG[nP][j].num=0;
      }
   }

   // QNORM ===============================================================
   if (p->Verbose) fprintf(stderr,"[1st-Step]");

   for (i=From; i< To; i++) { // Qnorm for each datafile: STEP 1
        LoadFile(fList, i, dataIn);
        if (p->Verbose) { fprintf(stderr,"."); fflush(stderr);}
        fprintf(stderr,"P8. tid=%d procesa exp i=%d\n",tid,i); fflush(stderr);

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
        fprintf(stderr,"P8.1tid=%d END exp i=%d\n",tid,i); fflush(stderr);

   } // end "for i" Qnorm for each datafile: STEP 1

  } // close Subsection [2b]

 #pragma omp section 
    { // open Subsection [2c]

         fprintf(stderr,"Pw. tid=%d Msj entre secciones\n",tid); fflush(stderr);
         fprintf(stderr,"Pz. tid=%d EN 2da SECC\n",tid); fflush(stderr);

    } // close Subsection [2c]


 } // close set of Sections [1]
} // close General parallel section [0]

   // HERE only one thread
   // Row average     [use col=tid=0] ------------------------------

     fprintf(stderr,"Px tid=%d ALONE totaliza\n",tid); fflush(stderr);
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
      fprintf(stderr, "P11 tid=%d checkVal=%lf\n",tid,checkVal); fflush(stderr);

      // Now copy the global averge into the local ones
      // *distribute* the global array
      for (i=0;i<nG;i++) {
        for (j=0;j<nP;j++) {
           AvG[j][i].Av =AvG[nP][i].Av;
           AvG[j][i].num=AvG[nP][i].num;
        }
      }


#ifdef DEBUG
     fprintf(stderr, "Row Average------------\n");
     for (j=0;j<nG;j++) fprintf (stderr,"%f (%d) ", AvG[j].Av,AvG[j].num); 
     fprintf(stderr,"\n");
#endif


   // Finally produce the ORDERED output file [STEP 2]-------------------------
   fprintf(stderr,"Here only tid=%d (must be Zero)\n",tid); fflush(stderr);

   if (!p->MemIndex) { 
      fclose(fI);
      if ((fI=fopen("~tmp","rb"))==NULL) terror("opening tmp-index for reading file");
   }

   // oputput file (by cols)
   if ((fOut=fopen(p->fOutName,"wb"))==NULL) {
      fprintf(stderr,"file : %s\n",p->fOutName); fflush(stderr);
       terror("opening OUTPUT file");
   }

   if ((dataOut=(double *)calloc(nG,sizeof(double)))==NULL) terror("memory for dataOut array");
   if (p->Verbose) { fprintf(stderr,"\n[2nd Step]"); fflush(stderr);}
   // NEWNEWNEWNEW
   if((dInde2=(int *)calloc(nG,sizeof(int)))==NULL) terror("memory for index2");

   fprintf(stderr,"P21 a grabar\n"); fflush(stderr);
 

   for (i=0;i<nE;i++) {
        if (p->Verbose) {fprintf(stderr,"."); fflush(stderr);}
        if (p->MemIndex) { // in memory - full
          for (j=0;j<nG;j++) { 
            dInde2[j]=mIndex[j][i];
          }
        } else {
          pos = fList[i].pos;
          fseek(fI, nG*pos*sizeof(int), SEEK_SET);
          fread(dInde2, sizeof(int), nG, fI);
        }

#ifdef DEBUG
        fprintf(stderr,"recovered Index (col=%d)\n",i);
        for (j=0;j<nG;j++) fprintf (stderr,"%d ", dIndex[j]); fprintf(stderr,"\n");
#endif

        // complete the output vector
        for (j=0;j<nG;j++) 
          dataOut[dInde2[j]]=AvG[tid][j].Av; // OJO

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
   fprintf(stderr,"P23 FIN FIN FINF IFN ....\n"); fflush(stderr);

   if (p->Traspose){
    if (p->Verbose) 
       { fprintf(stderr,"\nTransposing Disk->Disk"); fflush(stderr); }
    TransposeBin2Txt(p);
   }

   if (p->Verbose) fprintf(stderr,"Thread %d done.\n",tid);   

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

