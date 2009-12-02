/*******************************************************************************
  author       : Achim Gaedke
  filename     : ghmm/tests/coin_toss_test.c
  created      : DATE: 2001-04-25
  $Id: coin_toss_test.c 1521 2005-12-15 15:29:14Z cic99 $
*******************************************************************************/

#ifdef WIN32
#  include "win_config.h"
#endif

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdio.h>
#include <ghmm/matrix.h>
#include <ghmm/vector.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/model.h>
#include <ghmm/viterbi.h>
#include <ghmm/foba.h>

/*
  Simple model with one state and 2 symbols, like a coin toss
*/

int single_state_coin_toss()
{
  state single_state;
  model my_model;
  double symbols_single_state[2]={0.5,0.5};
  double trans_prob_single_state[1]={1.0};
  double trans_prob_single_state_rev[1]={1.0};
  int trans_id_single_state[1]={0};
  sequence_t* my_output;
  int silent_array[2] =  {0}; 

  my_model.model_type = 0;
  /* initialise this state */
  single_state.pi = 1.0;
  single_state.b=symbols_single_state;
  single_state.out_states=1;
  single_state.out_a=trans_prob_single_state;
  single_state.out_id=trans_id_single_state;
  single_state.in_states=1;
  single_state.in_id=trans_id_single_state;
  single_state.in_a=trans_prob_single_state_rev;
  single_state.fix=1;

  /* initialise model */
  my_model.N=1;
  my_model.M=2;
  my_model.s=&single_state;
  my_model.prior=-1;
  my_model.silent = silent_array;

  fprintf(stdout,"transition matrix:\n");
  model_A_print(stdout,&my_model,""," ","\n");
  fprintf(stdout,"observation symbol matrix:\n");
  model_B_print(stdout,&my_model,""," ","\n");

  my_output=model_generate_sequences(&my_model,0,10,10,100);
  sequence_print(stdout,my_output);

  sequence_free(&my_output);
  return 0;
}


/*
  another coin toss model
  flip between two states and output a different symbol
 */

int two_states_coin_toss()
{
  model my_model;
  state model_states[2];
  double symbols_head_state[2]={1.0,0.0};
  double trans_prob_head_state[2]={0.5,0.5};
  double trans_prob_head_state_rev[2]={0.5,0.5};
  int trans_id_head_state[2]={0,1};
  double symbols_tail_state[2]={0.0,1.0};
  double trans_prob_tail_state[2]={0.5,0.5};
  double trans_prob_tail_state_rev[2]={0.5,0.5};
  int trans_id_tail_state[2]={0,1};
  sequence_t *my_output;
  double log_p_viterbi, log_p_forward;
  double **forward_alpha;
  double forward_scale[10];
  int *viterbi_path;
  int i;
  /* flags indicating whether a state is silent */
  int silent_array[2] =  {0,0}; 

  my_model.model_type = 0;
  /* initialise head state */
  model_states[0].pi = 0.5;
  model_states[0].b=symbols_head_state;
  model_states[0].out_states=2;
  model_states[0].out_a=trans_prob_head_state;
  model_states[0].out_id=trans_id_head_state;
  model_states[0].in_states=2;
  model_states[0].in_id=trans_id_head_state;
  model_states[0].in_a=trans_prob_head_state_rev;
  model_states[0].fix=1;

  /* initialise tail state */
  model_states[1].pi = 0.5;
  model_states[1].b=symbols_tail_state;
  model_states[1].out_states=2;
  model_states[1].out_id=trans_id_tail_state;
  model_states[1].out_a=trans_prob_tail_state;
  model_states[1].in_states=2;
  model_states[1].in_id=trans_id_tail_state;
  model_states[1].in_a=trans_prob_tail_state_rev;
  model_states[1].fix=1;

  /* initialise model */
  my_model.N=2;
  my_model.M=2;
  my_model.s=model_states;
  my_model.prior=-1;
  my_model.silent = silent_array;
  
  fprintf(stdout,"transition matrix:\n");
  model_A_print(stdout,&my_model,""," ","\n");
  fprintf(stdout,"observation symbol matrix:\n");
  model_B_print(stdout,&my_model,""," ","\n");

  my_output=model_generate_sequences(&my_model,0,10,10,100);
  sequence_print(stdout,my_output);

  /* try viterbi algorithm in a clear situation */
  viterbi_path=viterbi(&my_model,
		       my_output->seq[0],
		       my_output->seq_len[0],
		       &log_p_viterbi);
  if (viterbi_path==NULL)
    {fprintf(stderr,"viterbi failed!"); return 1;}

  fprintf(stdout,"viterbi:\n");
  
  for(i=0;i<my_output->seq_len[0];i++){
    printf(" %d, ", viterbi_path[i]);
  }
  printf("\n");

  fprintf(stdout,
	  "log-p of this sequence (viterbi algorithm): %f\n",
	  log_p_viterbi);

  /* allocate matrix for forward algorithm */
  fprintf(stdout,"applying forward algorithm to the sequence...");
  forward_alpha=stat_matrix_d_alloc(10,2);
  if (forward_alpha==NULL)
    {
      fprintf(stderr,"\n could not alloc forward_alpha matrix\n");
      return 1;
    }

  /* run foba_forward */
  if (foba_forward(&my_model,
		   my_output->seq[0],
		   my_output->seq_len[0],
		   forward_alpha,
		   forward_scale,
		   &log_p_forward))
    {
      fprintf(stderr,"foba_logp failed!");
      stat_matrix_d_free(&forward_alpha);
      return 1;
    }

  /* alpha matrix */
  fprintf(stdout,"Done.\nalpha matrix from forward algorithm:\n");
  matrix_d_print(stdout,forward_alpha,10,2,""," ","\n");
  fprintf(stdout,"log-p of this sequence (forward algorithm): %f\n",log_p_forward);
  
  /* clean up */
  sequence_free(&my_output);
  free(viterbi_path);
  stat_matrix_d_free(&forward_alpha);
  return 0;
}

int main() {
  int result;

  /* Important! initialise rng  */
  ghmm_rng_init();

  if (single_state_coin_toss() || two_states_coin_toss())
    result = 1;
  else
    result = 0;
  
#ifdef WIN32
  printf("\nPress ENTER\n");
  fgetc(stdin);
#endif

  return result;
}
