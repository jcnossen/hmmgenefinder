/*******************************************************************************
  author       : Achim Gaedke
  filename     : ghmm/tests/coin_toss_test.c
  created      : DATE: 2001-04-25
  $Id: chmm.c 1195 2005-07-05 15:05:37Z schliep $
*******************************************************************************/

#ifdef WIN32
#  include "win_config.h"
#endif

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdio.h>
#include <unistd.h>
#include <ghmm/matrix.h>
#include <ghmm/rng.h>
#include <ghmm/smodel.h>

/*
  Simple model with one state and 2 symbols, like a coin toss
*/

int single_state_continuous()
{
  sstate single_state;
  smodel my_model;

  double trans_prob_single_state[]={1.0};
  double trans_prob_single_state_rev[]={1.0};
  double *trans_prob_single_state_array;
  double *trans_prob_single_state_rev_array;
  int trans_id_single_state[]={0};
  double c[]={1.0};
  double mue[]={0.0};
  double u[]={1.0};
  sequence_d_t* my_output;

  /* initialise transition array */
  trans_prob_single_state_array=trans_prob_single_state;
  trans_prob_single_state_rev_array=trans_prob_single_state_rev;

  /* initialise this state */
  single_state.pi = 1.0;
  single_state.out_states=1;
  single_state.out_a=&trans_prob_single_state_array;
  single_state.out_id=trans_id_single_state;
  single_state.in_states=1;
  single_state.in_id=trans_id_single_state;
  single_state.in_a=&trans_prob_single_state_rev_array;
  single_state.c=c;  /* weight of distribution */
  single_state.mue=mue; /* mean */
  single_state.u=u; /* variance */
  single_state.fix=0; /* training of output functions */

  /* initialise model */
  my_model.N=1; /* states */
  my_model.M=1; /* density functions per state */
  my_model.cos=1; /* class of states */
  my_model.density=0; /* normal distributions */
  my_model.prior=-1; /* a priori probability */
  my_model.s=&single_state; /* states array*/
 
#if 0
  /* print model */
  smodel_print(stdout,&my_model);
#endif

  /* generate sequences */
  my_output=smodel_generate_sequences(&my_model,
				      1,  /* random seed */
				      10, /* length of sequences */
				      10, /* sequences */
				      0,  /* label */
				      0  /* maximal sequence length 0: no limit*/
				      );
  /* print out sequences */
  sequence_d_print(stdout,    /* output file */
		   my_output, /* sequence */
		   0          /* do not truncate to integer*/
		   );


  /* reproduce the sequence by saving and rereading the model */
  {  
    FILE *my_file;
    sequence_d_t* new_output;
    char filename_buffer[]="/tmp/chmm_test.XXXXXX";
    int descriptor;
    int model_counter;
    smodel **model_array;

    descriptor=mkstemp(filename_buffer);
    /* write this model */
    my_file=fdopen(descriptor,"w+");
    fprintf(stdout,"printing model to file %s\n",filename_buffer);
    smodel_print(my_file,&my_model);
    (void)fseek(my_file, 0L, SEEK_SET);
    fclose(my_file);

    /* read this model */
    fprintf(stdout,"rereading model from file %s\n",filename_buffer);
    model_array=smodel_read(filename_buffer,&model_counter);

    /* generate sequences */
    fprintf(stdout,"generating sequences again\n");
    new_output=smodel_generate_sequences(model_array[0],
					 1,  /* random seed */
					 10, /* length of sequences */
					 10, /* sequences */
					 0,  /* label */
					 0   /* maximal sequence length 0: no limit*/
					 );

    sequence_d_print(stdout,     /* output file */
		     new_output, /* sequence */
		     0           /* do not truncate to integer*/
		     );
    /* free everything */
    close(descriptor);
    unlink(filename_buffer);
    sequence_d_free(&new_output);
    /*while(model_counter>0)
      {
	smodel_free(&(model_array[model_counter-1]));
	model_counter-=1;
	}*/
    smodel_free(model_array);
   
  }

  sequence_d_free(&my_output);
  return 0;
}

int main()
{
  /* Important! initialise rng  */
  ghmm_rng_init();

  return single_state_continuous();
}
