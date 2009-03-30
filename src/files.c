// exercise about how all threads can write in the same file 

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) 
{
int nthreads, tid;

/* Fork a team of threads giving them their own copies of variables */

#pragma omp parallel  private(nthreads, tid)
  {
  char txt[10];
  int i;
  FILE *f; // all threads have their own handle

  if ((f=fopen("borrar.txt","r+"))==NULL) // all threads open the file as Write/Read
     { printf("err abre file\n");exit(-1);}

  /* Obtain thread number */
  tid = omp_get_thread_num();

  for (i=0;i<10;i++) txt[i]=tid+'A';
  txt[9]='\n';

  nthreads = omp_get_num_threads();
  fseek(f,(tid*10),SEEK_SET);
  fprintf(f,"%s",txt);
  printf("1)thread %d ha escrito file (pos=%d)\n", tid,tid*10);


  fseek(f,((nthreads+tid)*10),SEEK_SET);
  fprintf(f,"%s",txt);
  printf("2)thread %d ha escrito file (pos=%d)\n", tid,(nthreads+tid)*10);
  fclose(f);

  }  /* All threads join master thread and disband */

}

