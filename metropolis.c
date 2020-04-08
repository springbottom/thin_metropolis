
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "metropolis.h"


int main(int argc, char *argv[]){

  if(argc != 7){
    printf("Not the correct number of arguments.");
  }

  //Importing command line arguments?
  int L        = atoi(argv[1]);
  int INIT     = atoi(argv[2]);
  double LOWER = strtod(argv[3],NULL);
  double UPPER = strtod(argv[4],NULL);
  int N        = atoi(argv[5]);
  int ITER     = atoi(argv[6]);

  //Setting up the file to print to?
  char file_str[200];
  sprintf(file_str,"%s_%s_%s_%s_%s_%s_%ld.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],time(NULL));
  fp = fopen(file_str,"w");

  //Actually doing things
  init();
  for (int i = 0; i < INIT){
    MCMC_step(1/LOWER);
  }
  for (double inverse_beta = LOWER; inverse_beta <= UPPER; inverse_beta += (UPPER-LOWER)/N.){
    fprintf(fp,"%f",inverse_beta);
    for (int i = 0; i < ITER){
      MCMC_step(1/inverse_beta);
    }
  }

  fprintf(fp,"test");
  fclose(fp);
  return 0;
}



double rand_double() {
  return rand()/(double)RAND_MAX;
}

int init_state(){
  for (int i = 0; i < 2; i++){
    for (int j = 0; j < L; j++){
      to_return[i][j][0] = cos(2*M_PI*rand_double());
      to_return[i][j][1] = sin(2*M_PI*rand_double());
    }
  }
  return 1;
}

int get_index(){
  index[0] = rand()%2;
  index[1] = rand()%2;

  n1[0] = (index[0]+1)%2; n1[1] = index[1];
  n2[0] = index[0]; n2[1] = (index[1]+1)%L;
  return 1;
}

int MCMC_step(double beta){
  get_index();

  double old_vector[2] = state[index[0]][index[1]];
  double vn1[2] = state[n1[0]][n1[1]];
  double vn2[2] = state[n2[0]][n2[1]];
  double neighbor[2] = {vn1[0]+vn2[0],vn1[1]+vn2[1]};

  double old_energy = old_vector[0]*neighbor[0]+old_vector[1]*neighbor[1];
  double new_angle = 2*M_PI*rand_double();
  double new_vector[2] = {cos(new_angle), sin(new_angle)};
  double new_energy = new_vector[0]*neighbor[0]+new_vector[1]*neighbor[1];

  if (rand_double() <= exp(beta*(old_energy-new_energy))){
    fprintf(fp,"%f,",new_energy);
    state[index[0]][index[1]] = new_vector;
    return 1;
  }
  fprintf(fp,"%f,",0.0);
  return 1;

  update()
}
