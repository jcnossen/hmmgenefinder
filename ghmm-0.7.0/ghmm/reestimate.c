/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/reestimate.c
*       Authors:  Bernhard Knab, Benjamin Georgi
*
*       Copyright (C) 1998-2004 Alexander Schliep
*       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
*       Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik,
*                               Berlin
*
*       Contact: schliep@ghmm.org
*
*       This library is free software; you can redistribute it and/or
*       modify it under the terms of the GNU Library General Public
*       License as published by the Free Software Foundation; either
*       version 2 of the License, or (at your option) any later version.
*
*       This library is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY; without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*       Library General Public License for more details.
*
*       You should have received a copy of the GNU Library General Public
*       License along with this library; if not, write to the Free
*       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*
*
*       This file is version $Revision: 1712 $
*                       from $Date: 2006-10-16 11:54:36 +0200 (Mon, 16 Oct 2006) $
*             last change by $Author: cic99 $.
*
*******************************************************************************/


#include "mes.h"
#include "mprintf.h"
#include "reestimate.h"
#include "matrix.h"
#include "model.h"
#include "foba.h"
#include "float.h"
#include "const.h"
#include "math.h"
#include "ghmm.h"
#include "kbestbasics.h"
#include <ghmm/internal.h>

typedef struct local_store_t {
  double *pi_num;
  double pi_denom;
  double **a_num;
  double *a_denom;
  double **b_num;
  double **b_denom;
} local_store_t;

static local_store_t *reestimate_alloc (const model * mo);
static int reestimate_free (local_store_t ** r, int N);
static int reestimate_init (local_store_t * r, const model * mo);
static int reestimate_setlambda (local_store_t * r, model * mo);
static int reestimate_one_step (model * mo, local_store_t * r,
                                int seq_number, int *seq_length, int **O,
                                double *log_p, double *seq_w);



static double nologSum(double* vec, int len) {
  int i;
  double retval=vec[0];
  
  for (i=1; i<len; i++)
    retval += vec[i];

  return retval;
}


/*----------------------------------------------------------------------------*/
static local_store_t *reestimate_alloc (const model * mo)
{
# define CUR_PROC "reestimate_alloc"
  int i;
  local_store_t *r = NULL;
  ARRAY_CALLOC (r, 1);

  ARRAY_CALLOC (r->pi_num, mo->N);
  ARRAY_CALLOC (r->a_num, mo->N);
  for (i = 0; i < mo->N; i++)
    ARRAY_CALLOC (r->a_num[i], mo->s[i].out_states);
  ARRAY_CALLOC (r->a_denom, mo->N);

  ARRAY_MALLOC (r->b_num, mo->N);
  for (i = 0; i < mo->N; i++)
    ARRAY_CALLOC (r->b_num[i], model_ipow(mo, mo->M, mo->s[i].order + 1));
  ARRAY_MALLOC (r->b_denom, mo->N);
  for (i = 0; i < mo->N; i++)
    ARRAY_CALLOC (r->b_denom[i], model_ipow (mo, mo->M, mo->s[i].order));

  return (r);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  reestimate_free (&r, mo->N);
  return (NULL);
# undef CUR_PROC
}                               /* reestimate_alloc */

/*----------------------------------------------------------------------------*/
static int reestimate_free (local_store_t ** r, int N)
{
# define CUR_PROC "reestimate_free"
  int i;
  mes_check_ptr (r, return (-1));
  if (!*r)
    return (0);
  m_free ((*r)->pi_num);

  if ((*r)->a_num){
    for (i = 0; i < N; i++){
      m_free ((*r)->a_num[i]);
    }  
  }  
  m_free ((*r)->a_num);
  m_free ((*r)->a_denom);

  if ((*r)->b_num){
    for (i = 0; i < N; i++){
      m_free ((*r)->b_num[i]);
    }
  }
  m_free ((*r)->b_num);
  if ((*r)->b_denom){
    for (i = 0; i < N; i++){
      m_free ((*r)->b_denom[i]);
    }  
  }
  m_free ((*r)->b_denom);

  m_free (*r);
  return (0);
# undef CUR_PROC
}                               /* reestimate_free */

/*----------------------------------------------------------------------------*/
static int reestimate_init (local_store_t * r, const model * mo)
{
# define CUR_PROC "reestimate_init"

  int i, j, m, size;

  r->pi_denom = 0.0;

  for (i = 0; i < mo->N; i++) {
    r->pi_num[i] = 0.0;
    r->a_denom[i] = 0.0;
    for (j = 0; j < mo->s[i].out_states; j++)
      r->a_num[i][j] = 0.0;

    size = model_ipow (mo, mo->M, mo->s[i].order);
    for (m = 0; m < size; m++)
      r->b_denom[i][m] = 0.0;
    size *= mo->M;
    for (m = 0; m < size; m++)
      r->b_num[i][m] = 0.0;
  }
  return (0);
# undef CUR_PROC
}                               /* reestimate_init */

/*----------------------------------------------------------------------------*/
int reestimate_alloc_matvek (double ***alpha, double ***beta, double **scale,
                             int T, int N)
{
# define CUR_PROC "reestimate_alloc_matvek"
  int res = -1;
  *alpha = stat_matrix_d_alloc (T, N);
  if (!(*alpha)) {
    mes_proc ();
    goto STOP;
  }
  *beta = stat_matrix_d_alloc (T, N);
  if (!(*beta)) {
    mes_proc ();
    goto STOP;
  }
  ARRAY_CALLOC (*scale, T);
  res = 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (res);
# undef CUR_PROC
}                               /* reestimate_alloc_matvek */

/*----------------------------------------------------------------------------*/
int reestimate_free_matvek (double **alpha, double **beta, double *scale, int T)
{
# define CUR_PROC "reestimate_free_matvek"
  stat_matrix_d_free (&alpha);
  stat_matrix_d_free (&beta);
  m_free (scale);
  return (0);
# undef CUR_PROC
}                               /* reestimate_free_matvek */

/*----------------------------------------------------------------------------*/
void reestimate_update_tie_groups (model * mo)
{
#define CUR_PROC "reestimate_update_tie_groups"
  int i, j, k;
  int bi_len;

  double *new_emissions;
  double nr = 0.0, non_silent_nr = 0.0;

  /* printf("** Start of reestimate_update_tie_groups **\n"); */

  /* do nothing if there are no tied emissions */
  if (!(mo->model_type & kTiedEmissions)) {
    printf ("No tied emissions in reestimate_update_tie_groups\n");
    return;
  }

  if (mo->model_type & kHigherOrderEmissions) {
    /*  printf("reestimate_update_tie_groups: Allocating for higher order states\n"); */
    ARRAY_MALLOC (new_emissions, model_ipow (mo, mo->M, mo->maxorder + 1));
  }
  else {
    /* printf("reestimate_update_tie_groups: No higher order states\n"); */
    ARRAY_MALLOC (new_emissions, mo->M);
  }

  for (i = 0; i < mo->N; i++) {
    bi_len = model_ipow (mo, mo->M, mo->s[i].order + 1);
    /* find tie group leaders */
    if (mo->tied_to[i] == i) {
      nr = 1.0;
      if (mo->silent[i] == 0) {
        non_silent_nr = 1.0;
      }
      else {
        non_silent_nr = 0.0;
      }

      /* printf("tie group leader %d found.\n",i); */

      /* initializing with tie group leader emissions */
      for (k = 0; k < bi_len; k++) {
        new_emissions[k] = mo->s[i].b[k];
      }

      /* finding tie group members */
      for (j = i + 1; j < mo->N; j++) {
        if (mo->tied_to[j] == i && mo->s[i].order == mo->s[j].order) {
          /* silent states have no contribution to the pooled emissions within a group */
          if (mo->silent[j] == 0) {
            nr += 1.0;
            non_silent_nr += 1.0;
            /* printf("  tie group member %d -> leader %d.\n",j,i); */
            /* summing up emissions in the tie group */
            for (k = 0; k < bi_len; k++) {
              new_emissions[k] += mo->s[j].b[k];
            }
          }
          /* updating silent flag */
          else {
            if (non_silent_nr > 0.0) {
              mo->silent[j] = 0;
            }
            nr += 1.0;
          }
        }
      }
      /*printf ("i = %d\n", i); */
      
        /* updating emissions */
      if (nr > 1.0) {
        for (j = i; j < mo->N; j++) {
          /* states within one tie group are required to have the same order */
          if (mo->tied_to[j] == i && mo->s[i].order == mo->s[j].order) {
            for (k = 0; k < bi_len; k++) {
              mo->s[j].b[k] = (new_emissions[k] / non_silent_nr);
              /* printf("s(%d)[%d] -> %f / %f = %f\n", j, k, new_emissions[k], nr,mo->s[j].b[k]);   */
            }
          }
        }
      }
    }
  }

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  m_free (new_emissions);
#undef CUR_PROC
}           /* reestimate_update_tie_groups */

/*----------------------------------------------------------------------------*/
static int reestimate_setlambda (local_store_t * r, model * mo)
{
# define CUR_PROC "reestimate_setlambda"
  int res = -1;
  int h, i, j, m, l, j_id, positive;
  double factor, p_i;
  int hist, col, size, reachable;
  char *str;

  mes_check_0 (r->pi_denom, goto STOP);

  for (i = 0; i < mo->N; i++) {
    reachable = 1;
    /* Pi */
    mo->s[i].pi = r->pi_num[i] / r->pi_denom;

    /* A */
    /* note: denom. might be 0; never reached state? */
    p_i = 0.0;
    if (r->a_denom[i] < EPS_PREC) {
      for (h = 0; h < mo->s[i].in_states; h++) {
        p_i += mo->s[i].in_a[h];
      }
      if (p_i == 0) {
        if (mo->s[i].in_states == 0)
          str = mprintf (NULL, 0,
                         " State %d can't be reached (no in_states)\n", i);
        else
          str =
            mprintf (NULL, 0, " State %d can't be reached (prob = 0)\n", i);
        mes_prot (str);
        m_free (str);
        reachable = 0;
      }
      factor = 0.0;
    }
    else
      factor = (1 / r->a_denom[i]);
    positive = 0;

    for (j = 0; j < mo->s[i].out_states; j++) {
      /* TEST: denom. < numerator */
      if ((r->a_denom[i] - r->a_num[i][j]) <= -EPS_PREC) {
        mes_prot (" numerator > denom!\n");
      }
      mo->s[i].out_a[j] = r->a_num[i][j] * factor;
      if (r->a_num[i][j] >= EPS_PREC)
        positive = 1;
      /* important: also update in_a  */
      l = 0;
      j_id = mo->s[i].out_id[j];
      while (l < mo->s[j_id].in_states)
        if (mo->s[j_id].in_id[l] == i)
          break;
        else
          l++;
      if (l == mo->s[j_id].in_states) {
        mes_prot ("no corresponding in_a to out_a found!\n");
        goto STOP;
      }
      mo->s[j_id].in_a[l] = mo->s[i].out_a[j];
    }

    /*if (!positive) {
       str = mprintf(NULL, 0, 
       "All numerator a[%d][j] == 0 (denom=%.4f, P(in)=%.4f)!\n", 
       i, r->a_denom[i], p_i);
       mes_prot(str);
       m_free(str);
       } */

    /* if fix, continue to next state */
    if (mo->s[i].fix)
      continue;

    /* B */
    size = model_ipow (mo, mo->M, mo->s[i].order);
    /* If all in_a's are zero, the state can't be reached.
       Set all b's to -1.0 */
    if (!reachable) {
      for (hist = 0; hist < size; hist++) {
	    col = hist * mo->M;
	    for (m = col; m < col + mo->M; m++) {
	      mo->s[i].b[m] = -1.0;
	    }
      }
    }
    else {
      for (hist = 0; hist < size; hist++) {
        /* If the denominator is very small, we have not seen many emissions
	   in this state with this history.
	   We are conservative and just skip them. */
	if (r->b_denom[i][hist] < EPS_PREC)
	  continue;
	else
	  factor = (1.0 / r->b_denom[i][hist]);

	positive = 0;
	/* TEST: denom. < numerator */
	col = hist * mo->M;
	for (m = col; m < col + mo->M; m++) {
	  if ((r->b_denom[i][hist] - r->b_num[i][m]) <= -EPS_PREC) {
	    str = mprintf (NULL, 0, "numerator b (%.4f) > denom (%.4f)!\n",
			   r->b_num[i][m], r->b_denom[i][hist]);
	    mes_prot (str);
	    m_free (str);
	  }
	  
	  mo->s[i].b[m] = r->b_num[i][m] * factor;
	  if (mo->s[i].b[m] >= EPS_PREC)
	    positive = 1;
	}

	if (!positive) {
	  str = mprintf (NULL, 0, "All numerator b[%d][%d-%d] == 0 (denom = %g)!\n",
			 i, col, col + mo->M, r->b_denom[i][hist]);
	  mes_prot (str);
	  m_free (str);
	}
      }                           /* for each history */
    }
  }                             /* for (i = 0 .. < mo->N)  */

  res = 0;
  if (mo->model_type & kTiedEmissions)
    reestimate_update_tie_groups (mo);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (res);
# undef CUR_PROC
}                               /* reestimate_setlambda */

/*----------------------------------------------------------------------------*/
static int reestimate_one_step (model * mo, local_store_t * r, int seq_number,
				int *seq_length, int **O, double *log_p,
				double *seq_w)
{
# define CUR_PROC "reestimate_one_step"
  int res = -1;
  int k, i, j, t, j_id, valid;
  int e_index;
  double **alpha = NULL;
  double **beta = NULL;
  double *scale = NULL;
  int T_k=0;
  double gamma;
  double log_p_k;

  /* first set maxorder to zero if model_type & kHigherOrderEmissions is FALSE 

     TODO XXX use model->maxorder only
     if model_type & kHigherOrderEmissions is TRUE */

  if (!(mo->model_type & kHigherOrderEmissions))
    mo->maxorder = 0;

  *log_p = 0.0;
  valid = 0;
  /* loop over all sequences */
  for (k = 0; k < seq_number; k++) {
    mo->emission_history = 0;
    T_k = seq_length[k];        /* current seq. length */

    /* initialization of  matrices and vector depends on T_k */
    if (reestimate_alloc_matvek (&alpha, &beta, &scale, T_k, mo->N) == -1) {
      mes_proc ();
      goto FREE;
    }
    if (foba_forward (mo, O[k], T_k, alpha, scale, &log_p_k) == -1) {
      mes_proc ();
      goto FREE;
    }

    if (log_p_k != +1) {        /* O[k] can be generated */
      *log_p += log_p_k;
      valid = 1;
      if (foba_backward (mo, O[k], T_k, beta, scale) == -1) {
        mes_proc ();
        goto FREE;
      }

      /* loop over all states */
      for (i = 0; i < mo->N; i++) {
        /* Pi */
        r->pi_num[i] += seq_w[k] * alpha[0][i] * beta[0][i];
        r->pi_denom += seq_w[k] * alpha[0][i] * beta[0][i];

        for (t = 0; t < T_k - 1; t++) {
          /* B */
          if (!mo->s[i].fix) {
            e_index = get_emission_index (mo, i, O[k][t], t);
            if (e_index != -1) {
              gamma = seq_w[k] * alpha[t][i] * beta[t][i];
              r->b_num[i][e_index] += gamma;
              r->b_denom[i][e_index / (mo->M)] += gamma;
            }
          }
          update_emission_history (mo, O[k][t]);

          /* A */
          r->a_denom[i] += (seq_w[k] * alpha[t][i] * beta[t][i]);
          for (j = 0; j < mo->s[i].out_states; j++) {
            j_id = mo->s[i].out_id[j];
            e_index = get_emission_index (mo, j_id, O[k][t + 1], t + 1);
            if (e_index != -1)
              r->a_num[i][j] += (seq_w[k] * alpha[t][i] * mo->s[i].out_a[j]
                                 * mo->s[j_id].b[e_index] * beta[t + 1][j_id]
                                 * (1.0 / scale[t + 1]));       /* c[t] = 1/scale[t] */
          }
        }
        /* B: last iteration for t==T_k-1 */
        t = T_k - 1;
        if (!mo->s[i].fix) {
          e_index = get_emission_index (mo, i, O[k][t], t);
          if (e_index != -1) {
            gamma = seq_w[k] * alpha[t][i] * beta[t][i];
            r->b_num[i][e_index] += gamma;
            r->b_denom[i][e_index / (mo->M)] += gamma;
          }
        }
      }
    }                           /* if (log_p_k != +1) */
    else
      printf ("O(%d) can't be built from model mo!\n", k);

    reestimate_free_matvek (alpha, beta, scale, T_k);

  }                             /* for (k = 0; k < seq_number; k++) */

  if (valid) {
    /* new parameter lambda: set directly in model */

    if (reestimate_setlambda (r, mo) == -1) {
      mes_proc ();
      goto STOP;
    }
    /* printf ("---- reestimate: after normalization ----\n"); */
    /*
       printf("Emission:\n");
       model_B_print(stdout, mo, "\t", " ", "\n");
     */
    /* only test: */
    if (model_check (mo) == -1) {
      mes_proc ();
      goto STOP;
    }
  }
  else {                        /* NO sequence can be built from model mo! */
    *log_p = +1;
  }

  res = 0;


STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
   return (res);
FREE:
   reestimate_free_matvek (alpha, beta, scale, T_k);
   return (res);
# undef CUR_PROC
}                               /* reestimate_one_step */


/*---------------------------------------------------------------------------*/
static int reestimate_one_step_lean (model * mo, local_store_t * r,
                                     int seq_number, int *seq_length,
                                     int **seqs, double *log_p, double *seq_w)
{
# define CUR_PROC "reestimate_one_step_lean"
  int res = -1;
  int i, t, k, e_index, len;
  int m, j, j_id, g, g_id, l, s, h;
  int * O;
  double c_t;
  
  double * alpha_last_col=NULL;
  double * alpha_curr_col=NULL;
  double * switching_tmp;
  double scale, old_scale, scalingf;
  double * summands=NULL;
  local_store_t * * curr_est=NULL, * * last_est=NULL, * * switch_lst;

  /* allocating memory for two columns of alpha matrix */
  ARRAY_CALLOC (alpha_last_col, mo->N);
  ARRAY_CALLOC (alpha_curr_col, mo->N);

  /* allocating 2*N local_store_t */
  ARRAY_CALLOC (last_est, mo->N);
  for (i=0; i<mo->N; i++)
    last_est[i] = reestimate_alloc(mo);
  ARRAY_CALLOC (curr_est, mo->N);
  for (i=0; i<mo->N; i++)
    curr_est[i] = reestimate_alloc(mo);


  /* temporary array to hold logarithmized summands
     for sums over probabilities */
  ARRAY_CALLOC (summands, m_max(mo->N,model_ipow(mo,mo->M,mo->maxorder+1))+1);


  for (k=0; k<seq_number; k++) {
    len = seq_length[k];
    O = seqs[k];
   
    foba_initforward(mo, alpha_last_col, O[0], &scale);
    if (scale < EPS_PREC) {
      /* means: first symbol can't be generated by hmm */
      *log_p = +1;
      goto STOP;
    }
    *log_p = - log(1/scale);
    
    for (t=1; t<len; t++) {
      old_scale = scale;
      scale = 0.0;
      update_emission_history(mo, O[t-1]);
    
      /* iterate over non-silent states */
      for (i=0; i<mo->N; i++) {
	/* printf("  akt_ state %d\n",i);*/
      
	e_index = get_emission_index(mo, i, O[t], t);
	if (e_index != -1){
	  alpha_curr_col[i] = foba_stepforward(&mo->s[i], alpha_last_col,
					       mo->s[i].b[e_index]);
	  scale += alpha_curr_col[i];
	}
	else
	  alpha_curr_col[i] = 0;
      }
    
      if (scale < EPS_PREC) {
	mes_prot("scale kleiner als eps\n");
	/* O-string  can't be generated by hmm */
	*log_p = +1.0;
	break;
      }
      c_t = 1/scale;
      for (i=0; i<mo->N; i++)
	alpha_curr_col[i] *= c_t;
    
      /*sum log(c[t]) scaling values to get  log( P(O|lambda) ) */
      *log_p -= log(c_t);

      scalingf = 1/old_scale;
      for (m=0; m<mo->N; m++) {
	for (i=0; i<mo->N; i++) {

	  /* computes estimates for the numerator of transition probabilities */
	  for (j=0; j<mo->s[i].out_states; j++) {
	    j_id = mo->s[i].out_id[j];

	    for (g=0; g<mo->s[j_id].in_states; g++) {
	      g_id = mo->s[j_id].out_id[g];
	      e_index = get_emission_index(mo, g_id, O[t], t);
	       /* scales all summands with the current */
	      summands[g] = last_est[m]->a_num[i][j] * mo->s[j_id].in_a[g]
		* mo->s[g_id].b[e_index] * scalingf;
	    }
	    if (j_id == m) {
	      e_index = get_emission_index(mo, j_id, O[t], t);
	      /* alpha is scaled. no other scaling necessary */
	      summands[g++] = alpha_last_col[i] * mo->s[i].out_a[j]
		* mo->s[j_id].b[e_index];
	    }
	    curr_est[m]->a_num[i][j] = nologSum(summands, g);
	  }
#ifdef calculate_denominators_extra
	  /* computes denominator of transition probabilities */	  
	  for (g=0; g<mo->s[m].in_states; g++) {
	    g_id = mo->s[m].in_id[g];
	    e_index = get_emission_index(mo, m, O[t], t);
	    /* scales all summands with the current factor */
	    summands[g] = last_est[m]->a_denom[i] * mo->s[m].in_a[g]
	      * mo->s[m].b[e_index] * scalingf;
	  }
	  if (i == m) {
	    for (l=0; l<mo->s[i].out_states; l++)
	      if (mo->s[i].out_id[l] == i)
		break;
	    if (l<mo->s[i].out_states) {
	      e_index = get_emission_index(mo, i, O[t], t);
	      /* alpha is scaled. no other scaling necessary */
	      summands[++g] = alpha_last_col[i] 
		* mo->s[i].out_a[l] * mo->s[m].b[e_index];
	    }
	  }
	  curr_est[m]->a_denom[i] = nologSum(summands, g);
#endif

	  /* computes estimates for the numerator of emmission probabilities*/
	  for (h=0; h<model_ipow(mo,mo->M, mo->s[i].order); h++)
	    for (s=h*mo->M; s<h*mo->M+mo->M; s++) {
	      for (g=0; g<mo->s[m].in_states; g++) {
		g_id = mo->s[m].out_id[g];
		e_index = get_emission_index(mo, g_id, O[t], t);
		/* scales all summands with the last scaling factor
		   of alpha */
		summands[g] = last_est[m]->b_num[i][s] * mo->s[m].in_a[g]
		  * mo->s[g_id].b[e_index] * scalingf;
	      }
	      curr_est[m]->b_num[i][s] = nologSum(summands, g);
	    }

	  e_index = get_emission_index(mo, m, O[t], t);
	  if (i == m) {
	    for (l=0; l<mo->s[m].out_states; l++)
	      if (mo->s[m].out_id[l] == m)
		break;
	    if (l<mo->s[m].out_states)
	      /* alpha is scaled. no other scaling necessary */
	      curr_est[m]->b_num[i][e_index] += alpha_last_col[i]
		* mo->s[m].out_a[l] * mo->s[m].b[e_index];
	  }
	  
	}
      }

      /* switching pointers of alpha_curr_col and alpha_last_col */
      switching_tmp  = alpha_last_col;
      alpha_last_col = alpha_curr_col;
      alpha_curr_col = switching_tmp;

      switch_lst = last_est;
      last_est   = curr_est;
      curr_est   = switch_lst;
    }

    /* filling the usual reestimate arrays by summing all states */
    for (m=0; m<mo->N; m++) {
      curr_est[m]->pi_denom = 0;
      for (i=0; i<mo->N; i++) {
	/* PI */
	/* XXX calculate the estimates for pi numerator */
	curr_est[m]->pi_num[i] = mo->s[i].pi;
	r->pi_num[i] += seq_w[k] * curr_est[m]->pi_num[i];
	curr_est[m]->pi_denom += curr_est[m]->pi_num[i];

	/* A */
	curr_est[m]->a_denom[i] = 0;
	for (j=0; j<mo->s[i].out_states; j++) {
	  r->a_num[i][j] += seq_w[k] * curr_est[m]->a_num[i][j];
	  curr_est[m]->a_denom[i] += curr_est[m]->a_num[i][j];
	}
	r->a_denom[i] += seq_w[k] * curr_est[m]->a_denom[i];
	
	/* B */
	for (h=0; h < model_ipow(mo,mo->M, mo->s[i].order); h++) {
	  curr_est[m]->b_denom[i][h] = 0;
	  for (s=h*mo->M; s<h*mo->M+mo->M; s++) {
	    r->b_num[i][s] += seq_w[k] * curr_est[m]->b_num[i][s];
	    curr_est[m]->b_denom[i][h] += curr_est[m]->b_num[i][s];
	  }
	  r->b_denom[i][h] += seq_w[k] * curr_est[m]->b_denom[i][h];
	}
      }
      /* PI */
      r->pi_denom += seq_w[k] * curr_est[m]->pi_denom;
    }
    
  }

  if (*log_p == 1.0)
    res = -1;
  else
    res = 0;

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  /* deallocation */
  if (last_est)
    for (i=0; i<mo->N; i++)
      reestimate_free(&(last_est[i]), mo->N);
  m_free(last_est);

  if (curr_est)
    for (i=0; i<mo->N; i++)
      reestimate_free(&(curr_est[i]), mo->N);
  m_free(curr_est);

  m_free(alpha_last_col);
  m_free(alpha_curr_col);
  m_free(summands);

  return res;  
#undef CUR_PROC
}


/*============================================================================*/
int reestimate_baum_welch (model * mo, sequence_t * sq)
{
# define CUR_PROC "reestimate_baum_welch"

  return reestimate_baum_welch_nstep (mo, sq, MAX_ITER_BW, EPS_ITER_BW);
# undef CUR_PROC
}                               /* reestimate_baum_welch */


/*============================================================================*/
int reestimate_baum_welch_nstep (model * mo, sequence_t * sq, int max_step,
                                 double likelihood_delta)
{
# define CUR_PROC "reestimate_baum_welch"
  int n, k, valid;
  double log_p, log_p_old, log_p_k, diff;
  local_store_t *r = NULL;
  int res = -1;
  char *str;

  /* local store for all iterations */
  r = reestimate_alloc (mo);
  if (!r) {
    mes_proc ();
    goto STOP;
  };

  log_p_old = -DBL_MAX;
  n = 1;

  /* main loop Baum-Welch-Alg. */
  while (n <= max_step) {
    printf("[i] Starting iteration %i / %i\n", n, max_step);
    if (1) {
      if (reestimate_one_step(mo, r, sq->seq_number, sq->seq_len, sq->seq,
			      &log_p, sq->seq_w) == -1) {
	str = mprintf (NULL, 0, "reestimate_one_step false (%d.step)\n", n);
	mes_prot (str);
	m_free (str);
	goto STOP;
      }
    }
    else {
      if (reestimate_one_step_lean(mo, r, sq->seq_number, sq->seq_len, sq->seq,
			      &log_p, sq->seq_w) == -1) {
	str = mprintf (NULL, 0, "reestimate_one_step false (%d.step)\n", n);
	mes_prot (str);
	m_free (str);
	goto STOP;
      }
    }

    /*if (n == 1)*/
    /*printf("%8.5f (-log_p input model)\n", -log_p); */
    /*else*/
    /*printf("%8.5f (-log_p)\n", -log_p); */

    if (log_p == +1) {
      printf
        ("Reestimate stopped: No sequence can be built from model mo!\n");
      break;
    }

    diff = log_p - log_p_old;
    /* error in convergence ? */
    if (diff < -EPS_PREC) {
      str = mprintf (NULL, 0, "No convergence: log P < log P-old! (n=%d)\n", n);
      mes_prot (str);
      m_free (str);
      goto STOP;
    }
    else if (log_p > EPS_PREC) {
      str = mprintf (NULL, 0, "No convergence: log P > 0! (n=%d)\n", n);
      mes_prot (str);
      m_free (str);
      goto STOP;
    }

    /* stop iterations? */
    if (diff < fabs ((double) likelihood_delta * log_p)) {
      printf("Convergence after %d steps\n", n); 

      break;
    }
    else {
      /* for next iteration */
      log_p_old = log_p;
      reestimate_init (r, mo);  /* sets all fields to zero */
      n++;
    }
  }                             /* while (n <= MAX_ITER) */

  /* log_p of reestimated model */
  log_p = 0.0;
  valid = 0;
  for (k = 0; k < sq->seq_number; k++) {
    if (foba_logp (mo, sq->seq[k], sq->seq_len[k], &log_p_k) == -1) {
      mes_proc ();
      goto STOP;
    }
    if (log_p_k != +1) {
      log_p += log_p_k;
      valid = 1;
    }
  }
  if (!valid)
    log_p = +1;
  /*printf("%8.5f (-log_p optimized model)\n", -log_p);*/

  /* check new parameter for plausibility */
  /* if (model_check(mo) == -1) { mes_proc(); goto STOP; } */
  res = 1;

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  reestimate_free (&r, mo->N);
  return res;
# undef CUR_PROC
}                               /* reestimate_baum_welch_nstep */





/**============================ Labeled HMMs ================================**/
/*----------------------------------------------------------------------------*/
static int reestimate_one_step_label (model * mo, local_store_t * r,
                                      int seq_number, int *seq_length,
                                      int **O, int **label, double *log_p,
                                      double *seq_w)
{
# define CUR_PROC "reestimate_one_step_label"

  int k, i, j, t, j_id;
  int e_index, T_k, valid = 0;
  double **alpha = NULL;
  double **beta = NULL;
  double *scale = NULL;
  double gamma;
  double log_p_k;
  char *str;

  /* first set maxorder to zero if model_type & kHigherOrderEmissions is FALSE 

     TODO XXX use model->maxorder only
     if model_type & kHigherOrderEmissions is TRUE */

  if (!(mo->model_type & kHigherOrderEmissions))
    mo->maxorder = 0;

  *log_p = 0.0;

  /* loop over all sequences */
  for (k = 0; k < seq_number; k++) {
    mo->emission_history = 0;
    T_k = seq_length[k];        /* current seq. length */

    /* initialization of  matrices and vector depends on T_k */
    if (reestimate_alloc_matvek (&alpha, &beta, &scale, T_k, mo->N) == -1) {
      mes_proc ();
      goto FREE;
    }

    if (foba_label_forward (mo, O[k], label[k], T_k, alpha, scale, &log_p_k)
        == -1) {
      mes_proc ();
      goto FREE;
    }

    if (log_p_k != +1) {        /* O[k] can be generated */
      *log_p += log_p_k;
      valid = 1;

      if (foba_label_backward (mo, O[k], label[k], T_k, beta, scale, &log_p_k)
          == -1) {
        mes_proc ();
        goto FREE;
      }

      /* loop over all states */
      for (i = 0; i < mo->N; i++) {
        /* Pi */
        r->pi_num[i] += seq_w[k] * alpha[0][i] * beta[0][i];
        r->pi_denom += seq_w[k] * alpha[0][i] * beta[0][i];

        for (t = 0; t < T_k - 1; t++) {
          /* B */
          if (!(mo->s[i].fix) && (mo->s[i].label == label[k][t])) {
            e_index = get_emission_index (mo, i, O[k][t], t);
            if (e_index != -1) {
              gamma = seq_w[k] * alpha[t][i] * beta[t][i];
              r->b_num[i][e_index] += gamma;
              r->b_denom[i][e_index / (mo->M)] += gamma;
            }
          }
          update_emission_history (mo, O[k][t]);

          /* A */
          r->a_denom[i] += seq_w[k] * alpha[t][i] * beta[t][i];
          for (j = 0; j < mo->s[i].out_states; j++) {
            j_id = mo->s[i].out_id[j];
            if (label[k][t + 1] != mo->s[j_id].label)
              continue;
            e_index = get_emission_index (mo, j_id, O[k][t + 1], t + 1);
            if (e_index != -1)
              r->a_num[i][j] += (seq_w[k] * alpha[t][i] * mo->s[i].out_a[j] *
                                 mo->s[j_id].b[e_index] * beta[t + 1][j_id] *
                                 (1.0 / scale[t + 1]));
          }
        }
        /* B: last iteration for t==T_k-1 */
        t = T_k - 1;
        if (!(mo->s[i].fix) && (mo->s[i].label == label[k][t])) {
          e_index = get_emission_index (mo, i, O[k][t], t);
          if (e_index != -1) {
            gamma = seq_w[k] * alpha[t][i] * beta[t][i];
            r->b_num[i][e_index] += gamma;
            r->b_denom[i][e_index / (mo->M)] += gamma;
          }
        }
      }
    }
    else {
      str = mprintf (NULL, 0, "warning: sequence %d can't be built from model\n", k);
      mes_prot (str);
      m_free (str);
    }

    reestimate_free_matvek (alpha, beta, scale, T_k);

  }                             /* for (k = 0; k < seq_number; k++) */

  if (valid) {
    /* new parameter lambda: set directly in model */
    if (reestimate_setlambda (r, mo) == -1) {
      mes_proc ();
      goto STOP;
    }
    /* printf ("---- reestimate: after normalization ----\n"); */
    /*     printf("Emission:\n"); */
    /*     model_B_print(stdout, mo, "\t", " ", "\n"); */
    /* only test: */
    /*    if (model_check(mo) == -1) { mes_proc(); goto STOP; } */
  }
  else {                        /* NO sequence can be built from model mo! */
    *log_p = +1;
  }

  return (0);

FREE:
  reestimate_free_matvek (alpha, beta, scale, T_k);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (-1);
# undef CUR_PROC
}                               /* reestimate_one_step_label */

/*============================================================================*/
int reestimate_baum_welch_label (model * mo, sequence_t * sq)
{
# define CUR_PROC "reestimate_baum_welch_label"

  return reestimate_baum_welch_nstep_label (mo, sq, MAX_ITER_BW, EPS_ITER_BW);
# undef CUR_PROC
}                               /* reestimate_baum_welch */

/*============================================================================*/
int reestimate_baum_welch_nstep_label (model * mo, sequence_t * sq,
                                       int max_step, double likelihood_delta)
{
# define CUR_PROC "reestimate_baum_welch_label"
  int n, k, valid;
  double log_p, log_p_old, log_p_k, diff;
  local_store_t *r = NULL;
  int res = -1;
  char *str;

  /* local store for all iterations */
  r = reestimate_alloc (mo);
  if (!r) {
    mes_proc ();
    goto STOP;
  };

  log_p_old = -DBL_MAX;
  n = 1;

  /* main loop Baum-Welch-Alg. */
  while (n <= max_step) {

    if (reestimate_one_step_label
        (mo, r, sq->seq_number, sq->seq_len, sq->seq, sq->state_labels,
         &log_p, sq->seq_w) == -1) {
      str = mprintf (NULL, 0, "reestimate_one_step_label false (%d.step)\n", n);
      mes_prot (str);
      m_free (str);
      goto STOP;
    }

    if (n == 1)
      printf ("%8.5f (-log_p input model)\n", -log_p);
    else
      printf ("%8.5f (-log_p)\n", -log_p);

    if (log_p == +1) {
      printf
        ("Reestimate stopped: No sequence can be built from model mo!\n");
      break;
    }

    diff = log_p - log_p_old;
    /* error in convergence ? */
    if (diff < -EPS_PREC) {
      str = mprintf (NULL, 0, "No convergence: log P < log P-old! (n = %d)\n", n);
      mes_prot (str);
      m_free (str);
      goto STOP;
    }
    else if (log_p > EPS_PREC) {
      str = mprintf (NULL, 0, "No convergence: log P > 0! (n = %d)\n", n);
      mes_prot (str);
      m_free (str);
      goto STOP;
    }

    /* stop iterations? */
    if (diff < fabs ((double) likelihood_delta * log_p)) {
      printf ("Convergence after %d steps\n", n);
      break;
    }
    else {
      /* for next iteration */
      log_p_old = log_p;
      reestimate_init (r, mo);  /* sets all fields to zero */
      n++;
    }
  }                             /* while (n <= MAX_ITER) */

  /* log_p of reestimated model */
  log_p = 0.0;
  valid = 0;
  for (k = 0; k < sq->seq_number; k++) {
    if (foba_label_logp (mo, sq->seq[k], sq->state_labels[k], sq->seq_len[k],
                         &log_p_k) == -1) {
      mes_proc ();
      goto STOP;
    }

    if (log_p_k != +1) {
      log_p += log_p_k;
      valid = 1;
    }
  }

  if (!valid)
    log_p = +1;
  printf ("%8.5f (-log_p optimized model)\n", -log_p);

  /* check new parameter for plausibility */
  /* if (model_check(mo) == -1) { mes_proc(); goto STOP; } */
  res = 0;

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  reestimate_free (&r, mo->N);
  return res;
# undef CUR_PROC
}                               /* reestimate_baum_welch_nstep_label */
