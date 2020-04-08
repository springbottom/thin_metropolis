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
  L        = atoi(argv[1]);
  INIT     = atoi(argv[2]);
  LOWER    = strtod(argv[3],NULL);
  UPPER    = strtod(argv[4],NULL);
  N        = atoi(argv[5]);
  ITER     = atoi(argv[6]);
  state    = malloc(L*4*sizeof(double));

  //Setting up the file to print to?
  char file_str[200];
  sprintf(file_str,"results/%s_%s_%s_%s_%s_%s_%ld.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],time(NULL));
  fp = fopen(file_str,"w");

  //Actually doing things
  init_state();
  for (double inverse_beta = LOWER; inverse_beta <= UPPER; inverse_beta += (UPPER-LOWER)/N){
    fprintf(fp,"%f\n",inverse_beta);
    for (int i = 0; i < INIT*2*L; i++){
      MCMC_step(1/inverse_beta);
    }
    energy = 0;
    for (int i = 0; i < ITER*2*L; i++){
      energy += MCMC_step(1/inverse_beta);
      //printf("%f,",energy);
      fprintf(fp,"%f,",energy);
    }
    fprintf(fp,"\n");
  }

  fclose(fp);
  free(state);
  return 0;
}

double rand_double() {
  return rand()/(double)RAND_MAX;
}

int get_index(i,j,k){
  return 4*i+2*j+k;
}

int init_state(){
  double this_theta = 2*M_PI*rand_double();
  for (int i = 0; i < L; i++){
    for (int j = 0; j < 2; j++){
      state[get_index(i,j,0)] = cos(this_theta);
      state[get_index(i,j,1)] = sin(this_theta);
    }
  }
  return 1;
}

/*
int get_random_index(){
  state_index[0] = rand()%L;
  state_index[1] = rand()%2;

  n1[0] = (state_index[0]+1)%L; n1[1] = state_index[1];
  n2[0] = state_index[0]; n2[1] = (state_index[1]+1)%2;
  n3[0] = (state_index[0]-1)%L; n3[1] = state_index[1];
  return 1;
}
*/

double MCMC_step(double beta){
  //get_random_index();

  int place = rand()%(2*L);//get_index(state_index[0],state_index[1],0);

  double *old_vector = &state[place];
  double *vn1 = &state[(place+4)%(2*L)];
  double *vn2 = &state[(place-4)%(2*L)];
  double *vn3 = &state[(place%2) == 2 ? place - 2 : place + 2];
  double neighbor[2] = {vn1[0]+vn2[0]+vn3[0],vn1[1]+vn2[1]+vn3[1]};

  double old_energy = old_vector[0]*neighbor[0]+old_vector[1]*neighbor[1];
  double new_angle = 2*M_PI*rand_double();
  double new_vector[2] = {cos(new_angle), sin(new_angle)};
  double new_energy = new_vector[0]*neighbor[0]+new_vector[1]*neighbor[1];

  if (rand_double() <= exp(beta*(old_energy-new_energy))){
    state[place]   = new_vector[0];
    state[place+1] = new_vector[1];
    return (new_energy - old_energy);
  }

  return 0;

}
