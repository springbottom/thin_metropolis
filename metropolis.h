

//Command Line ARGUMENTS
int L;         //the lattice size. (we work on (2,L) lattice)
int INIT;      //the number of initial iterations.
double LOWER;  //the lower bound for inverse beta
double UPPER;  //the upper bound for inverse beta
int N;         //the number of beta points to sample.
int ITER;      //the number of iterations per lattice point, per beta.

//Some variables
FILE *fp;
int state_index[2];
int n1[2];
int n2[2];
int n3[2];
double *state;
double energy;

//Function declarations
double H();
double rand_double();
int init_state();
int get_index();
double MCMC_step();
int positive_modulo(int i, int n);
