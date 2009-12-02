/*******************************************************************************
  author       : Janne Grunau
  filename     : ghmm/tests/label_higher_order_test.c
  created      : DATE: 2004-05-07
  $Id:
*******************************************************************************/

/* test_baumWelch
   generates a model to test c-functions with valgrind
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <ghmm/mes.h>
#include <ghmm/ghmm.h>
#include <ghmm/model.h>
#include <ghmm/viterbi.h>
#include <ghmm/reestimate.h>
#include <ghmm/gradescent.h>

#include <ghmm/matrix.h>
#include <ghmm/vector.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>

void generateModel(model *mo, int noStates) {
# define CUR_PROC "generateModel"

  state *states;
  int i, j;

  /* flags indicating whether a state is silent */
  int *silent_array;

  /*allocate memory for states and array of silent flags*/
  ARRAY_MALLOC(states, noStates);
  silent_array = (int*)malloc(sizeof(int)*noStates);

  /* initialize all states as none silent*/
  for (i=0; i < noStates; i++) {
    silent_array[i] = 0;
  }
 
  mo->N = noStates;
  mo->M = 4;
  mo->maxorder = noStates-1;
  mo->prior = -1;
  /* Model has Higher order Emissions and labeled states*/
  mo->model_type =  kLabeledStates;
  if (mo->maxorder>0)
    mo->model_type += kHigherOrderEmissions;
  /* kHigherOrderEmissions + kHasBackgroundDistributions*/

  /* allocate memory for pow look-up table and fill it */
  ARRAY_MALLOC(mo->pow_lookup, mo->maxorder+1)
  
  mo->pow_lookup[0] = 1;
  for (i=1; i<mo->maxorder+1; i++)
    mo->pow_lookup[i] =  mo->pow_lookup[i-1] * mo->M;

  /*initialize states*/
  for (i=0; i < mo->N; i++) {
    states[i].pi = (0==i ? 1.0:0.0);
    states[i].fix = 0;
    states[i].label = i%3;
    states[i].order = i%2;
    states[i].out_states = 2;
    states[i].in_states = 2;

    /* allocate memory for the a, the out- and incoming States and b array for higher emmission order states*/
    states[i].b = (double*)malloc(sizeof(double) * pow(mo->M, (states[i].order+1) ));
    states[i].out_id = (int*)malloc(sizeof(int)*states[i].out_states);
    states[i].in_id = (int*)malloc(sizeof(int)*states[i].in_states);
    states[i].out_a = (double*)malloc(sizeof(double)*states[i].out_states);
    states[i].in_a = (double*)malloc(sizeof(double)*states[i].in_states);

    for (j = 0; j < pow(mo->M,states[i].order+1); j++){
      states[i].b[j] = ( (0==(i+j)%mo->M) ? .6 : .4 / (mo->M-1));
    }

    if ((mo->N-1)==i) {
      states[i].out_id[0] = 0;
      states[i].out_id[1] = i;
    }
    else {
      states[i].out_id[0] = i;
      states[i].out_id[1] = i+1;
    }

    if (0==i) {
      states[i].in_id[0]  = i;
      states[i].in_id[1]  = mo->N-1;
    }
    else {
      states[i].in_id[1]  = i-1;
      states[i].in_id[0]  = i;
    }

    states[i].out_a[0] = 0.5;
    states[i].out_a[1] = 0.5;
    states[i].in_a[0]  = 0.5;
    states[i].in_a[1]  = 0.5;

#ifdef DEBUG
    printf("State %d goto    : %d, %d\n", i, states[i].out_id[0], states[i].out_id[1]);
    printf("State %d comefrom: %d, %d\n", i, states[i].in_id[0],  states[i].in_id[1]);
    printf("State %d goto    : %g, %g\n", i, states[i].out_a[0], states[i].out_a[1]);
    printf("State %d comefrom: %g, %g\n", i, states[i].in_a[0],  states[i].in_a[1]);
#endif
  }

  mo->s = states;
  mo->silent = silent_array;

#ifdef DEBUG
  for (i = 0; i < mo->N; i++) {
    printf("\n State %d:\n", i);
    for (j = 0; j < pow(mo->M,states[i].order+1); j++){
      printf("%g ",mo->s[i].b[j]);
    }
  }
#endif
  model_print(stdout, mo);

STOP:
  printf("\n");

# undef CUR_PROC
}


void testBaumwelch(){

  int i, error, tl,z,z1,z2;
  double log_p,first_prob1,first_prob2, first_prob;
  double *proba;
  int *path;
  int* real_path;
  int *path1;
  int* real_path1;
  int *path2;
  int* real_path2;
  model *mo = NULL;
  sequence_t *my_output, *your_output;
  int seqlen = 1000;
  tl = 150;

  mo = malloc(sizeof(model));
  if (mo==NULL) {fprintf(stderr,"Null Pointer in malloc(model).\n");}
  real_path = malloc(seqlen*sizeof(double));
  if(!real_path){ printf("real_path hat kein platz gekriegt\n");}
  real_path1 = malloc(seqlen*sizeof(double));
  if(!real_path1){ printf("real_path hat kein platz gekriegt\n");}
  real_path2 = malloc(seqlen*sizeof(double));
  if(!real_path2){ printf("real_path hat kein platz gekriegt\n");}
  /* generate a model with variable number of states*/
  generateModel(mo, 5);

  /*generate a random sequence*/
  my_output = model_label_generate_sequences(mo, 0, seqlen, 10, seqlen);
  for (i=0; i<seqlen; i++){
    printf("%d", my_output->state_labels[0][i]);
  }
  printf("\n");

  /*viterbi*/
  path = viterbi(mo, my_output->seq[0], my_output->seq_len[0], &first_prob);
  path1 = viterbi(mo, my_output->seq[1], my_output->seq_len[1], &first_prob1);
  path2 = viterbi(mo, my_output->seq[2], my_output->seq_len[2], &first_prob2);
  printf("\n viterbi-path\n");
  z=0;
  z1=0;
  z2=0;
  for (i=0; i<my_output->seq_len[0]; i++){
    if (path1[i] != -1) {
      real_path1[z1]=path1[i];
      z1++;
      printf("%d", path1[i]);
    }
    else printf("hallo");

    if (path2[i] != -1) {
      real_path2[z2]=path2[i];
      z2++;
      printf("%d", path2[i]);
    }
    else printf("hallo");

    if (path[i] != -1) {
      real_path[z]=path[i];
      z++;
      printf("%d", path[i]);
    }

  }
  printf("\n");
  printf("log-prob: %g\n",first_prob);
  my_output->state_labels[0]=real_path;
  my_output->state_labels[1]=real_path1;
  my_output->state_labels[2]=real_path2;

  for (i=0;i<seqlen;i++)
    printf("realpath[%i]=%i",i,real_path[i]);
  proba = malloc(sizeof(double)*tl);

  printf("No of Sequences = %d", my_output->seq_number);

  your_output = model_label_generate_sequences(mo, 0, seqlen, 1, seqlen);
  error = gradient_descent(&mo, your_output, .02, i);
  path = viterbi(mo, my_output->seq[0], my_output->seq_len[0], &log_p);
  free(path);

  /*reestimate_baum_welch_label(mo, my_output);*/
  /*reestimate_baum_welch(mo, my_output);*/

  /*reruns viterbi to check the training*/
  printf("run viterbi second\n");
  path = viterbi(mo, my_output->seq[0], my_output->seq_len[0], &log_p);
  for (i=0; i<(my_output->seq_len[0]*mo->N); i++){
    if (path[i] != -1) {printf("%d", path[i]);}
  }
  printf("\n");
  printf("log-prob: %g\n",log_p);


  /* freeing memory */
  model_free(&mo);
  free(path);
  /*printf("sequence_free success: %d\n", */sequence_free(&my_output)/*)*/;
  free(my_output);

}

int main(){

  ghmm_rng_init();
  testBaumwelch();

  return 0;
}
