/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/sgenerate.c
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
*       This file is version $Revision: 1276 $
*                       from $Date: 2005-08-24 17:34:29 +0200 (Wed, 24 Aug 2005) $
*             last change by $Author: cic99 $.
*
*******************************************************************************/

#ifdef WIN32
#  include "win_config.h"
#endif

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "mprintf.h"
#include "mes.h"
#include "const.h"
#include "sequence.h"
#include "smodel.h"
#include "sgenerate.h"
#include "sfoba.h"
#include "matrix.h"
#include "rng.h"
#include <ghmm/internal.h>


/*============================================================================*/

/* Extend given sequences on the basis of the model.
   Following modes are possible:
   mode = 0: Initial state: Viterbi, Extension: Viterbi-Path
        = 1: Initial state: Viterbi, Extension: all paths possible
	= 2: Initial state: probability(i) ~= alpha_t(i), 
	     Extension: Viterbi-Path
	= 3: Initial state: probability(i) ~= alpha_t(i),
	     Extension: all paths possible   
	     FRAGE: macht Extension Viterbi ueberhaupt Sinn???
*/



sequence_d_t *sgenerate_extensions (smodel * smo, sequence_d_t * sqd_short,
                                    int seed, int global_len,
                                    sgeneration_mode_t mode)
{
#define CUR_PROC "sgenerate_extensions"
  sequence_d_t *sq = NULL;
  int i, j, t, n, m, len = global_len, short_len, max_short_len = 0, up = 0;
#ifdef bausparkasse
  int tilgphase = 0;
#endif
  /* int *v_path = NULL; */
  double log_p, *initial_distribution, **alpha, *scale, p, sum;
  /* aicj */
  int class = -1;

  /* TEMP */
  if (mode == all_viterbi || mode == viterbi_viterbi || mode == viterbi_all) {
    mes_prot ("Error: mode not implemented yet\n");
    goto STOP;
  }

  if (len <= 0)
    /* no global length; model should have a final state */
    len = (int) MAX_SEQ_LEN;
  max_short_len = sequence_d_max_len (sqd_short);

  /*---------------alloc-------------------------------------------------*/
  sq = sequence_d_calloc (sqd_short->seq_number);
  if (!sq) {
    mes_proc ();
    goto STOP;
  }
  ARRAY_CALLOC (initial_distribution, smo->N);
  /* is needed in cfoba_forward() */
  alpha = matrix_d_alloc (max_short_len, smo->N);
  if (!alpha) {
    mes_proc ();
    goto STOP;
  }
  ARRAY_CALLOC (scale, max_short_len);
  ghmm_rng_init ();
  GHMM_RNG_SET (RNG, seed);

  /*---------------main loop over all seqs-------------------------------*/
  for (n = 0; n < sqd_short->seq_number; n++) {
    ARRAY_CALLOC (sq->seq[n], len);
    short_len = sqd_short->seq_len[n];
    if (len < short_len) {
      mes_prot ("Error: given sequence is too long\n");
      goto STOP;
    }
    sequence_d_copy (sq->seq[n], sqd_short->seq[n], short_len);
    sq->seq_label[n] = sqd_short->seq_label[n];

    /* Initial distribution */
    /* 1. Viterbi-state */
#if 0
    /* wieder aktivieren, wenn sviterbi realisiert */
    if (mode == viterbi_all || mode == viterbi_viterbi) {
      v_path = cviterbi (smo, sqd_short->seq[n], short_len, &log_p);
      if (v_path[short_len - 1] < 0 || v_path[short_len - 1] >= smo->N) {
        mes_prot ("Warning:Error: from viterbi()\n");
        sq->seq_len[n] = short_len;
        m_realloc (sq->seq[n], short_len);
        continue;
      }
      m_memset (initial_distribution, 0, smo->N);
      initial_distribution[v_path[short_len - 1]] = 1.0;        /* all other 0 */
      m_free (v_path);
    }
#endif

    /* 2. Initial Distribution ???
       Pi(i) = alpha_t(i)/P(O|lambda) */
    if (mode == all_all || mode == all_viterbi) {
      if (short_len > 0) {
        if (sfoba_forward (smo, sqd_short->seq[n], short_len, NULL /* ?? */ ,
                           alpha, scale, &log_p)) {
          mes_proc ();
          goto STOP;
        }
        sum = 0.0;
        for (i = 0; i < smo->N; i++) {
          /* alpha ist skaliert! */
          initial_distribution[i] = alpha[short_len - 1][i];
          sum += initial_distribution[i];
        }
        /* nicht ok.? auf eins skalieren? */
        for (i = 0; i < smo->N; i++)
          initial_distribution[i] /= sum;
      }
      else {
        for (i = 0; i < smo->N; i++)
          initial_distribution[i] = smo->s[i].pi;
      }
    }
    /* if short_len > 0:
       Initial state == final state from sqd_short; no output here
       else
       choose inittial state according to pi and do output
     */
    p = GHMM_RNG_UNIFORM (RNG);
    sum = 0.0;
    for (i = 0; i < smo->N; i++) {
      sum += initial_distribution[i];
      if (sum >= p)
        break;
    }
    /* error due to incorrect normalization ?? */
    if (i == smo->N) {
      i--;
      while (i > 0 && initial_distribution[i] == 0.0)
        i--;
    }
    t = 0;
    if (short_len == 0) {
      /* Output in state i */
      p = GHMM_RNG_UNIFORM (RNG);
      sum = 0.0;
      for (m = 0; m < smo->M; m++) {
        sum += smo->s[i].c[m];
        if (sum >= p)
          break;
      }
      /* error due to incorrect normalization ?? */
      if (m == smo->M) {
        m--;
        while (m > 0 && smo->s[i].c[m] == 0.0)
          m--;
      }
      sq->seq[n][t] = smodel_get_random_var (smo, i, m);

      if (smo->cos == 1) {
        class = 0;
      }
      else {
        if (!smo->class_change->get_class) {
          printf ("ERROR: get_class not initialized\n");
          goto STOP;
        }
        /*printf("1: cos = %d, k = %d, t = %d\n",smo->cos,smo->class_change->k,t);*/
        class = smo->class_change->get_class (smo, sq->seq[n], n, t);
      }


      t++;
    }
    /* generate completion for sequence */
    else {
      for (t = 0; t < short_len; t++)
        if (smo->cos == 1) {
          class = 0;
        }
        else {
          if (!smo->class_change->get_class) {
            printf ("ERROR: get_class not initialized\n");
            goto STOP;
          }
          /*printf("1: cos = %d, k = %d, t = %d\n",smo->cos,smo->class_change->k,t);*/
          class = smo->class_change->get_class (smo, sq->seq[n], n, t);
        }

      t = short_len;
    }
    while (t < len) {
      if (smo->s[i].out_states == 0)
        /* reached final state, exit while loop */
        break;
      sum = 0.0;
      for (j = 0; j < smo->s[i].out_states; j++) {
        sum += smo->s[i].out_a[class][j];
        if (sum >= p)
          break;
      }
      /* error due to incorrect normalization ?? */
      if (j == smo->s[i].out_states) {
        j--;
        while (j > 0 && smo->s[i].out_a[class][j] == 0.0)
          j--;
      }
      if (sum == 0.0) {
        /* Test: If an "empty" class, try the neighbour class;
           first, sweep down to zero, if still no success sweep up 
           to COs - 1. If still no success --> discard the sequence.
         */
        if (class > 0 && up == 0) {
          class--;
          continue;
        }
        else if (class < smo->cos - 1) {
          class++;
          up = 1;
          continue;
        }
        else {
          break;
        }
      }
      /* new state */
      i = smo->s[i].out_id[j];

      /* Output in state i */
      p = GHMM_RNG_UNIFORM (RNG);
      sum = 0.0;
      for (m = 0; m < smo->M; m++) {
        sum += smo->s[i].c[m];
        if (sum >= p)
          break;
      }
      if (m == smo->M) {
        m--;
        while (m > 0 && smo->s[i].c[m] == 0.0)
          m--;
      }
      /* random variable from density function */
      sq->seq[n][t] = smodel_get_random_var (smo, i, m);

      if (smo->cos == 1) {
        class = 0;
      }
      else {
        if (!smo->class_change->get_class) {
          printf ("ERROR: get_class not initialized\n");
          goto STOP;
        }
        /*printf("1: cos = %d, k = %d, t = %d\n",smo->cos,smo->class_change->k,t);*/
        class = smo->class_change->get_class (smo, sq->seq[n], n, t);
      }


      up = 0;
      t++;
    }                           /* while (t < len) */
    if (t < len)
      ARRAY_REALLOC (sq->seq[n], t);
    sq->seq_len[n] = t;

  }                             /* for n .. < seq_number */

  matrix_d_free (&alpha, max_short_len);
  m_free (scale);
  return sq;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  matrix_d_free (&alpha, max_short_len);
  sequence_d_free (&sq);
  return (NULL);
# undef CUR_PROC
}                               /* sgenerate_extensions */

/*============================================================================*/

double *sgenerate_single_ext (smodel * smo, double *O, const int len,
                              int *new_len, double **alpha,
                              sgeneration_mode_t mode)
{
# define CUR_PROC "sgenerate_single_ext"
  int i, j, m, t, class=0, up = 0;
  double *new_O = NULL, *scale = NULL, *initial_distribution = NULL;
  double log_p, sum, p;
  int k;
  /* TEMP */
  if (mode == all_viterbi || mode == viterbi_viterbi || mode == viterbi_all) {
    mes_prot ("Error: mode not implemented yet\n");
    goto STOP;
  }
  if (len <= 0) {
    mes_prot ("Error: sequence with zero or negativ length\n");
    goto STOP;
  }
  ARRAY_CALLOC (new_O, (int) MAX_SEQ_LEN);
  ARRAY_CALLOC (scale, len);
  ARRAY_CALLOC (initial_distribution, smo->N);
  sequence_d_copy (new_O, O, len);
  *new_len = len;
  /* Initial Distribution ???
     Pi(i) = alpha_t(i)/P(O|lambda) */
  if (mode == all_all || mode == all_viterbi) {
    if (sfoba_forward (smo, O, len, NULL /* ?? */ , alpha, scale, &log_p)) {
      mes_prot ("error from sfoba_forward, unable to extend\n");
      ARRAY_REALLOC (new_O, *new_len);
      return new_O;
    }
    sum = 0.0;


    for (i = 0; i < smo->N; i++) {
      /* alpha is scaled! */
      initial_distribution[i] = alpha[len - 1][i];
      sum += initial_distribution[i];
    }
    /* nicht ok.? scale to one? */
    for (i = 0; i < smo->N; i++) {
      initial_distribution[i] /= sum;
    }
  }

  p = GHMM_RNG_UNIFORM (RNG);
  sum = 0.0;
  for (i = 0; i < smo->N; i++) {
    sum += initial_distribution[i];
    if (sum >= p)
      break;
  }
  /* error due to incorrect normalization ?? */
  if (i == smo->N) {
    i--;
    while (i > 0 && initial_distribution[i] == 0.0)
      i--;
  }
  /* TEST */
  /* Already at the beginning in an end state? How is that possible? */
  if (smo->s[i].out_states == 0) {
    printf ("Beginn: Endzustand, State %d\n", i);
    for (k = 0; k < len; k++)
      printf ("%.2f ", O[k]);
    printf ("\n");
  }
  /* End Test */

  for (t = 0; t < len; t++)
    if (smo->cos == 1) {
      class = 0;
    }
    else {
      if (!smo->class_change->get_class) {
        printf ("ERROR: get_class not initialized\n");
        goto STOP;
      }
      /*printf("1: cos = %d, k = %d, t = %d\n",smo->cos,smo->class_change->k,t);*/
      class = smo->class_change->get_class (smo, O, 0, t);      /*XXX No sequence number */
    }

  t = len;
  while (t < (int) MAX_SEQ_LEN) {
    if (smo->s[i].out_states == 0)
      /* reached final state, exit while loop */
      break;
    p = GHMM_RNG_UNIFORM (RNG);
    sum = 0.0;
    for (j = 0; j < smo->s[i].out_states; j++) {
      sum += smo->s[i].out_a[class][j];
      if (sum >= p)
        break;
    }
    /* error due to incorrect normalization ?? */
    if (j == smo->s[i].out_states) {
      j--;
      while (j > 0 && smo->s[i].out_a[class][j] == 0.0)
        j--;
    }
    if (sum == 0.0) {
      /* Test: If an "empty" class, try the neighbour class;
         first, sweep down to zero, if still no success sweep up 
         to COs - 1. If still no success --> discard the sequence.
       */
      if (class > 0 && up == 0) {
        class--;
        continue;
      }
      else if (class < smo->cos - 1) {
        class++;
        up = 1;
        continue;
      }
      else {
        char *str =
          mprintf (NULL, 0, "unable to extend seq (all osc empty)\n");
        mes_prot (str);
        m_free (str);
        goto STOP;
      }
    }                           /* sum == 0 */

    /* new state */
    i = smo->s[i].out_id[j];

    if (smo->M == 1)
      m = 0;
    else {
      p = GHMM_RNG_UNIFORM (RNG);
      sum = 0.0;
      for (m = 0; m < smo->M; m++) {
        sum += smo->s[i].c[m];
        if (sum >= p)
          break;
      }
      if (m == smo->M) {
        m--;
        while (m > 0 && smo->s[i].c[m] == 0.0)
          m--;
      }
    }
    /* Output in state i, komp. m */
    /* random variable from density function */
    new_O[t] = smodel_get_random_var (smo, i, m);

    if (smo->cos == 1) {
      class = 0;
    }
    else {
      if (!smo->class_change->get_class) {
        printf ("ERROR: get_class not initialized\n");
        goto STOP;
      }
      /*printf("1: cos = %d, k = %d, t = %d\n",smo->cos,smo->class_change->k,t);*/
      class = smo->class_change->get_class (smo, new_O, 0, t);  /* XXX sequence number ? */
    }

    t++;
    up = 0;
  }                             /* while (t < MAX_SEQ_LEN) */
  if (t < (int) MAX_SEQ_LEN)
    ARRAY_REALLOC (new_O, t);

  *new_len = t;


  m_free (scale);
  m_free (initial_distribution);

  return new_O;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  m_free (new_O);
  m_free (scale);
  m_free (initial_distribution);
  return NULL;
# undef CUR_PROC
}                               /* sgenerate_single_ext */


/* generate a single next value bases on a trained model and on a seq und
   to length "len". Use the most prob. state given the seq as an initial state
   and determin the next state und the symbol with the RNG.

   If this function is called for one O successive "len", then it's possible
   to optimize some more.
*/

double sgenerate_next_value (smodel * smo, double *O, const int len)
{
# define CUR_PROC "sgenerate_next_value"
  double **alpha = NULL;
  double res = -1.0, sum, p;
  double *scale = NULL, log_p, max_val = -1000000;
  int i, j, m, init_state = -1;

  if (smo->cos > 1) {
    mes_prot ("sgenerate_next_value only for COS == 1\n");
    goto STOP;
  }

  alpha = matrix_d_alloc (len, smo->N);
  if (!alpha) {
    mes_proc ();
    goto STOP;
  }
  ARRAY_CALLOC (scale, len);
  if (sfoba_forward (smo, O, len, NULL /* ?? */ , alpha, scale, &log_p)) {
    mes_prot ("error from sfoba_forward\n");
    goto STOP;
  }

  /* find inititial state */
  sum = 0.0;
  for (i = 0; i < smo->N; i++)
    sum += alpha[len - 1][i];
  if (sum < 0.9 || sum > 1.1) {
    printf ("Error sum = %.4f (!= 1)\n", sum);
    goto STOP;
  }
  /* max state */
  for (i = 0; i < smo->N; i++) {
    if (alpha[len - 1][i] > max_val) {
      init_state = i;
      max_val = alpha[len - 1][i];
    }
  }

  /* random state */
  /*
     p = GHMM_RNG_UNIFORM(RNG);
     sum = 0.0;
     for (i = 0; i < smo->N; i++) {
     sum += alpha[len - 1][i];
     if (sum >= p)
     break;    
     }   
     if (i == smo->N) {
     i--;
     while (i > 0 && alpha[len - 1][i] == 0.0) i--;
     }  
     init_state = i;
   */


  if (init_state == -1 || smo->s[init_state].out_states == 0)
    goto STOP;

  p = GHMM_RNG_UNIFORM (RNG);
  sum = 0.0;
  for (j = 0; j < smo->s[init_state].out_states; j++) {
    sum += smo->s[init_state].out_a[0][j];
    if (sum >= p)
      break;
  }
  /* error due to incorrect normalization ?? */
  if (j == smo->s[init_state].out_states) {
    j--;
    while (j > 0 && smo->s[init_state].out_a[0][j] == 0.0)
      j--;
  }

  /* new state */
  i = smo->s[init_state].out_id[j];

  if (smo->M == 1)
    m = 0;
  else {
    p = GHMM_RNG_UNIFORM (RNG);
    sum = 0.0;
    for (m = 0; m < smo->M; m++) {
      sum += smo->s[i].c[m];
      if (sum >= p)
        break;
    }
    if (m == smo->M) {
      m--;
      while (m > 0 && smo->s[i].c[m] == 0.0)
        m--;
    }
  }
  /* Output in state i, komp. m */
  /* random variable from density function */
  res = smodel_get_random_var (smo, i, m);

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  matrix_d_free (&alpha, len);
  m_free (scale);
  return res;
# undef CUR_PROC
}
