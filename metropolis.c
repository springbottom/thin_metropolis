#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "metropolis.h"


int main(int argc, char *argv[]){

  //Importing command line arguments?
  if(argc != 7){
    printf("Not the correct number of arguments. Loading some defaults.\n");
    L        = 8;
    INIT     = 10000;
    LOWER    = 0.2;
    UPPER    = 1.5;
    N        = 10;
    ITER     = 5;
  }
  else{
    L        = atoi(argv[1]);
    INIT     = atoi(argv[2]);
    LOWER    = strtod(argv[3],NULL);
    UPPER    = strtod(argv[4],NULL);
    N        = atoi(argv[5]);
    ITER     = atoi(argv[6]);
  }

  //Setting up the file to print to?
  char file_str[200];
  sprintf(file_str,"results/%i_%i_%.1f_%.1f_%i_%i_%ld.txt",L,INIT,LOWER,UPPER,N,ITER,time(NULL));
  fp = fopen(file_str,"w");

  //Actually doing things
  state    = (double*)malloc(L*4*sizeof(double));
  init_state();
  for (double inverse_beta = LOWER; inverse_beta <= UPPER; inverse_beta += (UPPER-LOWER)/N){
    fprintf(fp,"%f\n",inverse_beta);
    for (int i = 0; i < INIT*2*L; i++){
      MCMC_step(1/inverse_beta);
    }
    energy = H();
    //printf("%f\n",energy);
    for (int i = 0; i < ITER*2*L; i++){
      energy += MCMC_step(1/inverse_beta);
      //printf("%f,%f\n",energy,H());
      fprintf(fp,"%f,",energy);
    }
    printf("%f,%f\n",energy,H());
    fprintf(fp,"\n");
  }

  fclose(fp);
  free(state);
  return 0;
}

inline int positive_modulo(int i, int n) {
    return (i % n + n) % n;
}

double rand_double() {
  return rand()/(double)RAND_MAX;
}

int get_index(i,j,k){
  return 4*i+2*j+k;
}

int init_state(){
  for (int i = 0; i < L; i++){
    for (int j = 0; j < 2; j++){
      double this_theta = 2*M_PI*rand_double();
      state[get_index(i,j,0)] = cos(this_theta);
      state[get_index(i,j,1)] = sin(this_theta);
    }
  }
  return 1;
}

//Mainly for debugging purposes
double H(){
  double to_return = 0;
  int ind1;
  for (int i = 0; i < L; i++){
    ind1 = get_index(i,0,0);
    to_return += state[ind1]*state[(ind1+4)%(4*L)] + state[(ind1+1)%(4*L)]*state[(ind1+5)%(4*L)];
    to_return += state[ind1]*state[(ind1+2)%(4*L)] + state[(ind1+1)%(4*L)]*state[(ind1+3)%(4*L)];
    to_return += state[(ind1+2)%(4*L)]*state[(ind1+6)%(4*L)] + state[(ind1+3)%(4*L)]*state[(ind1+7)%(4*L)];
    //printf("%f\n",to_return);
  }
  return -to_return;
}

double MCMC_step(double beta){

  int place = 2*(rand()%(2*L));
  //printf("%d\n",place);

  double *old_vector = &state[place];
  double *vn1 = &state[(place+4)%(4*L)];
  double *vn2 = &state[positive_modulo(place-4,4*L)]; //this only matters on the left edge...
  double *vn3 = &state[(place%4) == 2 ? (place - 2) : (place + 2)]; //actually I dont need the mods here
  double neighbor[2] = {vn1[0]+vn2[0]+vn3[0],vn1[1]+vn2[1]+vn3[1]};

  double old_energy = -old_vector[0]*neighbor[0]-old_vector[1]*neighbor[1];
  double new_angle = 2*M_PI*rand_double();
  double new_vector[2] = {cos(new_angle), sin(new_angle)};
  double new_energy = -new_vector[0]*neighbor[0]-new_vector[1]*neighbor[1];

  if (rand_double() <= exp(beta*(old_energy-new_energy))){
    state[place]   = new_vector[0];
    state[place+1] = new_vector[1];
    return (new_energy - old_energy);
  }

  return 0;

}
