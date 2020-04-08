
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double rand_double() {
  return rand()/(double)RAND_MAX;
}

double new_vector[2];

double * update(FILE* fp, double beta, double old_vector[2], double neighbor[2]) {
  double old_energy = old_vector[0]*neighbor[0]+old_vector[1]*neighbor[1];
  double new_angle = 2*M_PI*rand_double();
  new_vector[0] = cos(new_angle); new_vector[1] = sin(new_angle);
  double new_energy = new_vector[0]*neighbor[0]+new_vector[1]*neighbor[1];

  if (rand_double() <= exp(beta*(old_energy-new_energy))){
    fprintf(fp,"%f",new_energy);
    return new_vector;
  }
  fprintf(fp,"%f",0.0);
  return old_vector;
}


int main(int argc, char *argv[]){
  /*
  ARGUMENTS:
    -[1] : L, the lattice size. (we work on (2,N) lattice)
    -[2] : INIT, the number of initial iterations.
    -[3] : LOWER, the lower bound for beta
    -[4] : UPPER, the upper bound for beta
    -[5] : N, the number of beta points to sample.
    -[6] : ITER, the number of iterations per lattice point, per beta.
  */

  if(argc != 7){
    printf("Not the correct number of arguments.");
  }

  //Importing command line arguments?
  int L = atoi(argv[1]);
  int INIT = atoi(argv[2]);
  char*pEND;
  double LOWER = strtod(argv[3],&pEND);
  double UPPER = strtod(argv[4],&pEND);
  int N = atoi(argv[5]);
  int ITER = atoi(argv[6]);

  //Setting up the file to print to?
  FILE *fp;
  char file_str[200];
  sprintf(file_str,"%s_%s_%s_%s_%s_%s_%ld.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],time(NULL));
  fp = fopen(file_str,"w");

  fprintf(fp,"test");
  fclose(fp);
  return 0;
}
