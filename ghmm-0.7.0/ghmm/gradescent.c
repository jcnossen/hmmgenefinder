/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/gradescent.c
*       Authors:  Janne Grunau, Alexander Riemer
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
*       This file is version $Revision: 1293 $
*                       from $Date: 2005-09-02 19:32:09 +0200 (Fri, 02 Sep 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mes.h"
#include "matrix.h"
#include "model.h"
#include "foba.h"
#include "reestimate.h"
#include "gradescent.h"
#include "ghmm.h"
#include <ghmm/internal.h>

static double compute_performance (model * mo, sequence_t * sq);

/*----------------------------------------------------------------------------*/
static void gradient_descent_gfree (double **matrix_b, double *matrix_a,
                             double *matrix_pi, int N)
{
#define CUR_PROC "gradient_descent_gfree"

  int i;

  if (matrix_b)
    for (i = 0; i < N; i++) {
      m_free (matrix_b[i]);
    }
  m_free (matrix_b);

  m_free (matrix_a);
  m_free (matrix_pi);

#undef CUR_PROC
}


/*----------------------------------------------------------------------------*/
/** allocates memory for m and n matrices: */
static int gradient_descent_galloc (double ***matrix_b, double **matrix_a,
                             double **matrix_pi, model * mo)
{
#define CUR_PROC "gradient_descent_galloc"

  int i;

  /* first allocate memory for matrix_b */
  ARRAY_MALLOC (*matrix_b, mo->N);
  for (i = 0; i < mo->N; i++)
    ARRAY_CALLOC ((*matrix_b)[i], model_ipow (mo, mo->M, mo->s[i].order + 1));

  /* matrix_a(i,j) = matrix_a[i*mo->N+j] */
  ARRAY_CALLOC (*matrix_a, mo->N * mo->N);

  /* allocate memory for matrix_pi */
  ARRAY_CALLOC (*matrix_pi, mo->N);

  return 0;

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  gradient_descent_gfree (*matrix_b, *matrix_a, *matrix_pi, mo->N);
  return -1;

#undef CUR_PROC
}


/*----------------------------------------------------------------------------*/
/**
   computes matrices of n and m variables (expected values for how often a
   certain parameter from A or B is used)
   computes Baum-Welch variables implicit 
   @return                 nothing
   @param mo:              pointer to a model
   @param alpha:           matrix of forward variables
   @param backward:        matrix of backward variables
   @param scale:           scaling vector from forward-backward-algorithm
   @param seq:             sequence in internal representation
   @param seq_len:         length of sequence
   @param matrix_b:        matrix for parameters from B (n_b or m_b)
   @param matrix_a:        matrix for parameters from A (n_a or m_a)
   @param vec_pi:          vector for parameters in PI (n_pi or m_pi)
 */
int gradescent_compute_expectations (model * mo, double **alpha,
                                     double **beta, double *scale, int *seq,
                                     int seq_len, double **matrix_b,
                                     double *matrix_a, double *vec_pi)
{
#define CUR_PROC "gradescent_compute_expectations"

  int h, i, j, t;

  int size, j_id, e_index;

  double gamma, xi;
  double foba_sum;

  /* initialise matrices with zeros */
  for (i = 0; i < mo->N; i++) {
    for (j = 0; j < mo->N; j++)
      matrix_a[i * mo->N + j] = 0;
    size = model_ipow (mo, mo->M, mo->s[i].order + 1);
    for (h = 0; h < size; h++)
      matrix_b[i][h] = 0;
  }

  for (t = 0; t < seq_len; t++) {

    /* sum products of forward and backward variables over all states: */
    foba_sum = 0.0;
    for (i = 0; i < mo->N; i++)
      foba_sum += alpha[t][i] * beta[t][i];
    if (EPS_PREC > fabs (foba_sum)) {
      printf
        ("gradescent_compute_expect: foba_sum (%g) smaller as EPS_PREC (%g). t = %d.\n",
         foba_sum, EPS_PREC, t);
      return -1;
    }

    for (i = 0; i < mo->N; i++) {

      /* compute gamma implicit */
      gamma = alpha[t][i] * beta[t][i] / foba_sum;

      /* n_pi is easiest: n_pi(i) = gamma(0,i) */
      if (0 == t)
        vec_pi[i] = gamma;

      /* n_b(i,c) = sum[t, hist(t)=c | gamma(t,i)] / sum[t | gamma(t,i)] */
      e_index = get_emission_index (mo, i, seq[t], t);
      if (-1 != e_index)
        matrix_b[i][e_index] += gamma;
    }

    /* updating history, xi needs the right e_index for the next state */
    update_emission_history (mo, seq[t]);

    for (i = 0; i < mo->N; i++) {
      /* n_a(i,j) = sum[t=0..T-2 | xi(t,i,j)] / sum[t=0..T-2 | gamma(t,i)] */
      /* compute xi only till the state before the last */
      for (j = 0; (j < mo->s[i].out_states) && (t < seq_len - 1); j++) {
        j_id = mo->s[i].out_id[j];

        /* compute xi implicit */
	xi = 0;
        e_index = get_emission_index (mo, j_id, seq[t + 1], t + 1);
        if (e_index != -1)
          xi = alpha[t][i] * beta[t + 1][j_id] * mo->s[i].out_a[j]
            * mo->s[j_id].b[e_index] / (scale[t + 1] * foba_sum);

        matrix_a[i * mo->N + j_id] += xi;
      }
    }
  }

  return 0;
#undef CUR_PROC
}


/*----------------------------------------------------------------------------*/
static double compute_performance (model * mo, sequence_t * sq)
{
#define CUR_PROC "compute_performance"

  int k, seq_len, success = 1;
  /* log P[O | lambda, labeling] as computed by the forward algorithm */
  double log_p;
  /* sum over log P (calculated by forward_label) for all sequences
     used to compute the performance of the training */
  double log_p_sum = 0.0;

  /* loop over all sequences */
  for (k = 0; k < sq->seq_number && success; k++) {
    success = 0;
    seq_len = sq->seq_len[k];

    if (-1 !=
        foba_label_logp (mo, sq->seq[k], sq->state_labels[k], seq_len,
                         &log_p)) {
      success = 1;
      log_p_sum += log_p;

      if (-1 != foba_logp (mo, sq->seq[k], seq_len, &log_p))
        log_p_sum -= log_p;
      else {
        printf ("foba_forward error in sequence %d, length: %d\n", k,
                seq_len);
        success = 0;
      }
    }
    else
      printf ("foba_label_forward error in sequence %d, length: %d\n", k,
              seq_len);
  }

  /* return log_p_sum in success or +1.0 a probality of 0.0 on error */
  if (success)
    return log_p_sum;
  else
    return 1.0;
#undef CUR_PROC
}


/*----------------------------------------------------------------------------*/
/**
   Trains the model with a set of annotated sequences using gradient descent.
   Model must not have silent states. (one iteration)
   @return            0/-1 success/error
   @param mo:         pointer to a model
   @param sq:         struct of annotated sequences
   @param eta:        training parameter for gradient descent
 */
static int gradient_descent_onestep (model * mo, sequence_t * sq, double eta)
{
#define CUR_PROC "gradient_descent_onestep"

  int g, h, i, j, k;
  int seq_len, j_id, hist, size;

  /* expected usage of model parameters (A and B entries) */
  double *m_pi, *n_pi, *m_a, *n_a;
  double **m_b, **n_b;

  /* forward & backward variables w/ their scaling vector */
  double **alpha, **beta, *scale;

  /* variables to store values associated with the algorithm */
  double pi_sum, a_row_sum, b_block_sum, gradient;

  /* log P[O | lambda, labeling] as computed by the forward algorithm */
  double log_p;

  /* allocate memory for the parameters used for reestimation */
  if (-1 == gradient_descent_galloc (&m_b, &m_a, &m_pi, mo))
    return -1;
  if (-1 == gradient_descent_galloc (&n_b, &n_a, &n_pi, mo))
    return -1;

  /* loop over all sequences */
  for (k = 0; k < sq->seq_number; k++) {
    seq_len = sq->seq_len[k];

    if (-1 == reestimate_alloc_matvek (&alpha, &beta, &scale, seq_len, mo->N))
      continue;

    /* calculate forward and backward variables without labels: */
    if (-1 == foba_forward (mo, sq->seq[k], seq_len, alpha, scale, &log_p)) {
      printf ("forward error!\n");
      goto FREE;
    }

    if (-1 == foba_backward (mo, sq->seq[k], seq_len, beta, scale)) {
      printf ("backward error!\n");
      goto FREE;
    }

    /* compute n matrices (no labels): */
    if (-1 ==
        gradescent_compute_expectations (mo, alpha, beta, scale, sq->seq[k],
                                         seq_len, m_b, m_a, m_pi))
      printf ("Error in sequence %d, length %d (no labels)\n", k, seq_len);

    /* calculate forward and backward variables with labels: */
    if (-1 ==
        foba_label_forward (mo, sq->seq[k], sq->state_labels[k], seq_len,
                            alpha, scale, &log_p)) {
      printf ("forward labels error!\n");
      goto FREE;
    }
    if (-1 ==
        foba_label_backward (mo, sq->seq[k], sq->state_labels[k], seq_len,
                             beta, scale, &log_p)) {
      printf ("backward labels error!\n");
      goto FREE;
    }

    /* compute m matrices (labels): */
    if (-1 ==
        gradescent_compute_expectations (mo, alpha, beta, scale, sq->seq[k],
                                         seq_len, m_b, m_a, m_pi))
      printf ("Error in sequence %d, length %d (with labels)\n", k, seq_len);


    /* reestimate model parameters: */
    /* PI */
    pi_sum = 0;
    /*  update */
    for (i = 0; i < mo->N; i++) {

      if (mo->s[i].pi > 0.0) {
        gradient = eta * (m_pi[i] - n_pi[i]);
        if (mo->s[i].pi + gradient > EPS_PREC)
          mo->s[i].pi += gradient;
        else
          mo->s[i].pi = EPS_PREC;
      }

      /* sum over new PI vector */
      pi_sum += mo->s[i].pi;
    }
    if (pi_sum < EPS_PREC) {
      /* never get here */
      fprintf (stderr, "Training ruined the model. You lose.\n");
      k = sq->seq_number;
      goto FREE;
    }
    /*  normalise */
    for (i = 0; i < mo->N; i++)
      mo->s[i].pi /= pi_sum;

    /* A */
    for (i = 0; i < mo->N; i++) {
      a_row_sum = 0;
      /* update */
      for (j = 0; j < mo->s[i].out_states; j++) {
        j_id = mo->s[i].out_id[j];

        gradient =
          eta * (m_a[i * mo->N + j_id] - n_a[i * mo->N + j_id]) / (seq_len -
                                                                   1);
        if (mo->s[i].out_a[j] + gradient > EPS_PREC)
          mo->s[i].out_a[j] += gradient;
        else
          mo->s[i].out_a[j] = EPS_PREC;

        /* sum over rows of new A matrix */
        a_row_sum += mo->s[i].out_a[j];
      }

      if (a_row_sum < EPS_PREC) {
        /* never get here */
        fprintf (stderr, "Training ruined the model. You lose.\n");
        k = sq->seq_number;
        goto FREE;
      }
      /* normalise */
      for (j = 0; j < mo->s[i].out_states; j++) {
        mo->s[i].out_a[j] /= a_row_sum;

        /* mirror out_a to corresponding in_a */
        j_id = mo->s[i].out_id[j];
        for (g = 0; g < mo->s[j_id].in_states; g++)
          if (i == mo->s[j_id].in_id[g]) {
            mo->s[j_id].in_a[g] = mo->s[i].out_a[j];
            break;
          }
      }
    }

    /* B */
    for (i = 0; i < mo->N; i++) {

      /* don't update fix states */
      if (mo->s[i].fix)
        continue;

      /* update */
      size = model_ipow (mo, mo->M, mo->s[i].order);
      for (h = 0; h < size; h++) {
        b_block_sum = 0;
        for (g = 0; g < mo->M; g++) {
          hist = h * mo->M + g;
          gradient = eta * (m_b[i][hist] - n_b[i][hist]) / seq_len;
          /* printf("gradient[%d][%d] = %g, m_b = %g, n_b = %g\n"
             , i, hist, gradient, m_b[i][hist], n_b[i][hist]); */
          if (gradient + mo->s[i].b[hist] > EPS_PREC)
            mo->s[i].b[hist] += gradient;
          else
            mo->s[i].b[hist] = EPS_PREC;

          /* sum over M-length blocks of new B matrix */
          b_block_sum += mo->s[i].b[hist];
        }
        if (b_block_sum < EPS_PREC) {
          /* never get here */
          fprintf (stderr, "Training ruined the model. You lose.\n");
          k = sq->seq_number;
          goto FREE;
        }
        /* normalise */
        for (g = 0; g < mo->M; g++) {
          hist = h * mo->M + g;
          mo->s[i].b[hist] /= b_block_sum;
        }
      }
    }

    /* restore "tied_to" property */
    if (mo->model_type & kTiedEmissions)
      reestimate_update_tie_groups (mo);

  FREE:
    reestimate_free_matvek (alpha, beta, scale, seq_len);
  }

  gradient_descent_gfree (m_b, m_a, m_pi, mo->N);
  gradient_descent_gfree (n_b, n_a, n_pi, mo->N);

  return 0;

#undef CUR_PROC
}


/*----------------------------------------------------------------------------*/
/**
   Trains the model with a set of annotated sequences till convergence using
   gradient descent.
   Model must not have silent states. (checked in Python wrapper)
   @return            0/-1 success/error
   @param mo:         pointer to a model
   @param sq:         struct of annotated sequences
   @param eta:        intial parameter eta (learning rate)
   @param no_steps    number of training steps
 */
int gradient_descent (model ** mo, sequence_t * sq, double eta, int no_steps)
{
#define CUR_PROC "gradient_descent_nstep"

  int runs = 0;
  double cur_perf, last_perf;
  model *last;

  last = (model *) model_copy (*mo);
  last_perf = compute_performance (last, sq);

  while (eta > EPS_PREC && runs < no_steps) {
    runs++;
    if (-1 == gradient_descent_onestep (*mo, sq, eta)) {
      model_free (&last);
      return -1;
    }
    cur_perf = compute_performance (*mo, sq);

    if (last_perf < cur_perf) {
      /* if model is degenerated, lower eta and try again */
      if (cur_perf > 0.0) {
        printf ("current performance = %g\n", cur_perf);
        model_free (mo);
        *mo = (model *) model_copy (last);
        eta *= .5;
      }
      else {
        /* Improvement insignificant, assume convergence */
        if (fabs (last_perf - cur_perf) < cur_perf * (-1e-8)) {
          model_free (&last);
          printf ("convergence after %d steps.\n", runs);
          return 0;
        }

        if (runs < 175 || 0 == runs % 50)
          printf ("Performance: %g\t improvement: %g\t step %d\n", cur_perf,
                  cur_perf - last_perf, runs);

        /* significant improvement, next iteration */
        model_free (&last);
        last = (model *) model_copy (*mo);
        last_perf = cur_perf;
        eta *= 1.07;
      }
    }
    /* no improvement */
    else {

      if (runs < 175 || 0 == runs % 50)
        printf ("Performance: %g\t !IMPROVEMENT: %g\t step %d\n", cur_perf,
                cur_perf - last_perf, runs);

      /* try another training step */
      runs++;
      eta *= .85;
      if (-1 == gradient_descent_onestep (*mo, sq, eta)) {
        model_free (&last);
        return -1;
      }
      cur_perf = compute_performance (*mo, sq);
      printf ("Performance: %g\t ?Improvement: %g\t step %d\n", cur_perf,
              cur_perf - last_perf, runs);
      /* improvement, save and proceed with next iteration */
      if (last_perf < cur_perf && cur_perf < 0.0) {
        model_free (&last);
        last = (model *) model_copy (*mo);
        last_perf = cur_perf;
      }
      /* still no improvement, revert to saved model */
      else {
        runs--;
        model_free (mo);
        *mo = (model *) model_copy (last);
        eta *= .9;
      }
    }
  }

  model_free (&last);
  return 0;

#undef CUR_PROC
}
