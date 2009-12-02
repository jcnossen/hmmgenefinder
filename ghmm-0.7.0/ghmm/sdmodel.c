/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/sdmodel.c
*       Authors:  Wasinee Rungsarityotin, Utz Pape, Andrea Weisse
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
*       This file is version $Revision: 1264 $
*                       from $Date: 2005-08-10 18:33:54 +0200 (Wed, 10 Aug 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/

#ifdef WIN32
#  include "win_config.h"
#endif

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#undef NDEBUG
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "sdmodel.h"
#include "model.h"
#include "matrix.h"
#include "sequence.h"
#include "const.h"
#include "rng.h"
#include "sdfoba.h"
#include "mes.h"
#include "mprintf.h"
#include "string.h"
#include "ghmm.h"
#include <ghmm/internal.h>

#define  __EPS 10e-6

#ifdef sdmodelSTATIC
/*----------------------------------------------------------------------------*/
static int sdmodel_state_alloc (sdstate * state, int M, int in_states,
                                int out_states, int cos)
{
#define CUR_PROC "sdmodel_state_alloc"
  int res = -1;
  ARRAY_CALLOC (state->b, M);

  if (out_states > 0) {
    ARRAY_CALLOC (state->out_id, out_states);
    state->out_a = matrix_d_alloc (cos, out_states);
    if (!state->out_a) {
      mes_proc ();
      goto STOP;
    }
  }
  if (in_states > 0) {
    ARRAY_CALLOC (state->in_id, in_states);
    state->in_a = matrix_d_alloc (cos, in_states);
    if (!state->in_a) {
      mes_proc ();
      goto STOP;
    }
  }

  res = 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (res);
#undef CUR_PROC
}                               /* model_state_alloc */
#endif

/*----------------------------------------------------------------------------*/
double sdmodel_likelihood (sdmodel * mo, sequence_t * sq)
{
# define CUR_PROC "sdmodel_likelihood"
  double log_p_i, log_p;
  int found, i;

  found = 0;
  log_p = 0.0;
  for (i = 0; i < sq->seq_number; i++) {
    sdfoba_logp (mo, sq->seq[i], sq->seq_len[i], &log_p_i);

    if (log_p_i != +1) {
      log_p += log_p_i;
      found = 1;
    }
    /*    else {*/
    /*  char *str =*/
    /*   mprintf(NULL, 0, "sequence[%d] can't be build.\n", i);*/
    /*  mes_prot(str);*/
    /*}*/
  }
  if (!found)
    log_p = +1.0;
  return (log_p);
#undef CUR_PROC
}                               /*sdmodel_likelihood */

#ifdef sdmodelSTATIC
/*----------------------------------------------------------------------------*/
static int sdmodel_copy_vectors (sdmodel * mo, int index, double ***a_matrix,
                                 double **b_matrix, double *pi, int *fix)
{
#define CUR_PROC "sdmodel_copy_vectors"
  int i, c, cnt_out = 0, cnt_in = 0;

  mo->s[index].pi = pi[index];
  mo->s[index].fix = fix[index];
  for (i = 0; i < mo->M; i++)
    mo->s[index].b[i] = b_matrix[index][i];
  for (c = 0; c < mo->cos; c++)
    for (i = 0; i < mo->N; i++) {
      if (a_matrix[c][index][i]) {      /* Transitions to a following state possible */
        if (cnt_out >= mo->s[index].out_states) {
          mes_proc ();
          return (-1);
        }
        mo->s[index].out_id[cnt_out] = i;
        mo->s[index].out_a[c][cnt_out] = a_matrix[c][index][i];
        cnt_out++;
      }
      if (a_matrix[i][index]) { /* Transitions to a previous state possible */
        if (cnt_in >= mo->s[index].in_states) {
          mes_proc ();
          return (-1);
        }
        mo->s[index].in_id[cnt_in] = i;
        mo->s[index].in_a[c][cnt_in] = a_matrix[c][i][index];
        cnt_in++;
      }
    }
  return (0);
#undef CUR_PROC
}                               /* model_alloc_vectors */
#endif

/*============================================================================*/

/* Old prototyp:

model **model_read(char *filename, int *mo_number, int **seq,
			 const int *seq_len, int seq_number) { */



/*============================================================================*/

int sdmodel_free (sdmodel ** mo)
{
#define CUR_PROC "sdmodel_free"
  sdstate *my_state;
  int i;
  mes_check_ptr (mo, return (-1));
  if (!*mo)
    return (0);
  for (i = 0; i < (*mo)->N; i++) {
    my_state = &((*mo)->s[i]);
    if (my_state->b)
      m_free (my_state->b);
    if (my_state->out_id != NULL)
      m_free (my_state->out_id);
    if (my_state->in_id != NULL)
      m_free (my_state->in_id);
    /*
       if (my_state->out_a)
       m_free(my_state->out_a);
       if (my_state->in_a)
       m_free(my_state->in_a); */
    if (my_state->out_a)
      matrix_d_free (&((*mo)->s[i].out_a), (*mo)->cos);
    if (my_state->in_a)
      matrix_d_free (&((*mo)->s[i].in_a), (*mo)->cos);
    my_state->pi = 0;
    my_state->b = NULL;
    my_state->out_id = NULL;
    my_state->in_id = NULL;
    my_state->out_a = NULL;
    my_state->in_a = NULL;
    my_state->out_states = 0;
    my_state->in_states = 0;
    my_state->fix = 0;
    m_free (my_state->label);
  }
  m_free ((*mo)->s);
  m_free (*mo);
  fprintf (stderr, "Free sdmodel\n");
  return (0);
#undef CUR_PROC
}                               /* model_free */


/*============================================================================*/
sdmodel *sdmodel_copy (const sdmodel * mo)
{
# define CUR_PROC "sdmodel_copy"
  int i, j, k, nachf, vorg, m;
  sdmodel *m2 = NULL;
  ARRAY_CALLOC (m2, 1);
  ARRAY_CALLOC (m2->s, mo->N);
  for (i = 0; i < mo->N; i++) {
    nachf = mo->s[i].out_states;
    vorg = mo->s[i].in_states;
    ARRAY_CALLOC (m2->s[i].out_id, nachf);
    m2->s[i].out_a = matrix_d_alloc (mo->cos, nachf);
    ARRAY_CALLOC (m2->s[i].in_id, vorg);
    m2->s[i].in_a = matrix_d_alloc (mo->cos, vorg);

    ARRAY_CALLOC (m2->s[i].b, mo->M);
    /* Copy the values */

    for (j = 0; j < nachf; j++) {
      for (k = 0; k < mo->cos; k++) {
        m2->s[i].out_a[k][j] = mo->s[i].out_a[k][j];
      }
      m2->s[i].out_id[j] = mo->s[i].out_id[j];
    }
    for (j = 0; j < vorg; j++) {
      for (k = 0; k < mo->cos; k++) {
        m2->s[i].in_a[k][j] = mo->s[i].in_a[k][j];
      }
      m2->s[i].in_id[j] = mo->s[i].in_id[j];
    }

    for (m = 0; m < mo->M; m++) {
      m2->s[i].b[m] = mo->s[i].b[m];
    }
    m2->s[i].pi = mo->s[i].pi;
    m2->s[i].out_states = nachf;
    m2->s[i].in_states = vorg;
    m2->s[i].label = (char *) malloc (strlen (mo->s[i].label) + 1);
    strcpy (m2->s[i].label, mo->s[i].label);
    m2->s[i].countme = mo->s[i].countme;
  }
  m2->N = mo->N;
  m2->M = mo->M;
  m2->prior = mo->prior;
  m2->cos = mo->cos;
  m2->model_type = mo->model_type;
  if (mo->model_type == kSilentStates) {
    assert (mo->silent != NULL);
    ARRAY_CALLOC (m2->silent, mo->N);
    for (i = 0; i < mo->N; i++) {
      m2->silent[i] = mo->silent[i];
    }
    if (mo->topo_order_length > 0) {
      ARRAY_CALLOC (m2->topo_order, mo->topo_order_length);
      for (i = 0; i < mo->topo_order_length; i++) {
        m2->topo_order[i] = mo->topo_order[i];
      }
    }
  }
  if (mo->get_class) {
    m2->get_class = mo->get_class;
  }
  return (m2);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sdmodel_free (&m2);
  return (NULL);
# undef CUR_PROC
}                               /* model_copy */


/*----------------------------------------------------------------------------*/
void sdmodel_topo_ordering (sdmodel * mo)
{
#define CUR_PROC "sdmodel_topo_ordering"
  fprintf (stderr, "sdmodel_topo_ordering will be implemented using DFS.\n");
#undef CUR_PROC
}

#ifdef sdmodelSTATIC
/*============================================================================*/
static sequence_t *__sdmodel_generate_sequences (sdmodel * mo, int seed,
                                                 int global_len,
                                                 long seq_number, int Tmax)
{
# define CUR_PROC "sdmodel_generate_sequences"

  /* An end state is characterized by not having an output probabiliy. */
  unsigned long tm;             /* Time seed */
  sequence_t *sq = NULL;
  int state, n, i, j, m, reject_os, reject_tmax, badseq, class;
  double p, sum, osum = 0.0;
  int len = global_len, up = 0, stillbadseq = 0, reject_os_tmp = 0;
  double dummy = 0.0;

  sq = sequence_calloc (seq_number);
  if (!sq) {
    mes_proc ();
    goto STOP;
  }
  if (len <= 0)
    /* A specific length of the sequences isn't given. As a model should have
       an end state, the constant MAX_SEQ_LEN is used. */
    len = (int) MAX_SEQ_LEN;

  if (seed > 0) {
    GHMM_RNG_SET (RNG, seed);
  }
  else {
    tm = rand ();
    GHMM_RNG_SET (RNG, tm);
    fprintf (stderr, "# using rng '%s' seed=%ld\n", GHMM_RNG_NAME (RNG), tm);
  }

  n = 0;
  reject_os = reject_tmax = 0;

  while (n < seq_number) {
    /* Test: A new seed for each sequence */
    /*   ghmm_rng_timeseed(RNG); */
    stillbadseq = badseq = 0;
    ARRAY_CALLOC (sq->seq[n], len);

    /* Get a random initial state i */
    p = GHMM_RNG_UNIFORM (RNG);
    sum = 0.0;
    for (i = 0; i < mo->N; i++) {
      sum += mo->s[i].pi;
      if (sum >= p)
        break;
    }

    /* Get a random initial output m */
    p = GHMM_RNG_UNIFORM (RNG);
    sum = 0.0;
    for (m = 0; m < mo->M; m++) {
      sum += mo->s[i].b[m];
      if (sum >= p)
        break;
    }
    sq->seq[n][0] = m;
    state = 1;

    /* The first symbol chooses the start class */
    class = mo->get_class (sq->seq[n], state);
    /*class = sequence_d_class(&dummy, 0, &osum); */ /*  dummy function */
    while (state < len) {

      /* Get a new state */
      p = GHMM_RNG_UNIFORM (RNG);
      sum = 0.0;
      for (j = 0; j < mo->s[i].out_states; j++) {
        sum += mo->s[i].out_a[class][j];
        if (sum >= p)
          break;
      }

      if (sum == 0.0) {
        if (mo->s[i].out_states > 0) {
          /* Repudiate the sequence, if all smo->s[i].out_a[class][.] == 0,
             that is, class "class" isn't used in the original data:
             go out of the while-loop, n should not be counted. */
          /* printf("Zustand %d, class %d, len %d out_states %d \n", i, class,
             state, smo->s[i].out_states); */
          badseq = 1;
          /* break; */

          /* Try: If the class is "empty", try the neighbour class;
             first, sweep down to zero; if still no success, sweep up to
             COS - 1. If still no success --> Repudiate the sequence. */
          if (class > 0 && up == 0) {
            class--;
            continue;
          }
          else if (class < mo->cos - 1) {
            class++;
            up = 1;
            continue;
          }
          else {
            stillbadseq = 1;
            break;
          }
        }
        else
          /* An end state is reached, get out of the while-loop */
          break;
      }
      i = mo->s[i].out_id[j];

      /* Get a random output m from state i */
      p = GHMM_RNG_UNIFORM (RNG);
      sum = 0.0;
      for (m = 0; m < mo->M; m++) {
        sum += mo->s[i].b[m];
        if (sum >= p)
          break;
      }

      sq->seq[n][state] = m;

      /* Decide the class for the next step */
      class = mo->get_class (sq->seq[n], state);
      /*class = sequence_d_class(&dummy, state, &osum); */ /* dummy */
      up = 0;
      state++;
    }                           /* while (state < len) , global_len depends on the data */

    if (badseq) {
      reject_os_tmp++;
    }

    if (stillbadseq) {
      reject_os++;
      m_free (sq->seq[n]);
      /*      printf("cl %d, s %d, %d\n", class, i, n); */
    }
    else if (state > Tmax) {
      reject_tmax++;
      m_free (sq->seq[n]);
    }
    else {
      if (state < len)
        ARRAY_REALLOC (sq->seq[n], state);
      sq->seq_len[n] = state;
      /* sq->seq_label[n] = label; */
      /* vector_d_print(stdout, sq->seq[n], sq->seq_len[n]," "," ",""); */
      n++;
    }
    /*    printf("reject_os %d, reject_tmax %d\n", reject_os, reject_tmax); */
    if (reject_os > 10000) {
      mes_prot ("Reached max. no. of rejections\n");
      break;
    }
    if (!(n % 1000))
      printf ("%d Seqs. generated\n", n);
  }                             /* n-loop */


  if (reject_os > 0)
    printf ("%d sequences rejected (os)!\n", reject_os);
  if (reject_os_tmp > 0)
    printf ("%d sequences changed class\n", reject_os_tmp - reject_os);
  if (reject_tmax > 0)
    printf ("%d sequences rejected (Tmax, %d)!\n", reject_tmax, Tmax);

  return (sq);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_free (&sq);
  return (NULL);
# undef CUR_PROC
}                               /* data */
#endif

int sdmodel_initSilentStates (sdmodel * mo)
{
#define CUR_PROC "sdmodel_initSilentStates"
  int nSilentStates = 0;
  int i, m;
  double sum;
  int *__silents;

  ARRAY_CALLOC (__silents, mo->N);

  for (i = 0; i < mo->N; i++) {
    for (m = 0, sum = 0; m < mo->M; m++) {
      sum += mo->s[i].b[m];
    }
    if (sum < __EPS) {
      __silents[i] = 1;
      nSilentStates++;
    }
    else
      __silents[i] = 0;
  }

  if (nSilentStates) {
    mo->model_type = kSilentStates;
    mo->silent = __silents;
  }
  else {
    mo->model_type = kNotSpecified;
    mo->silent = NULL;
    m_free (__silents);
  }
  return (nSilentStates);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (0);
#undef CUR_PROC
}

/*============================================================================*/
/*
 * Before calling generate sequence, the global random number generator must
 * be initialized by calling ghmm_rng_init().
 *
 */
/*returns the extended sequence struct with state matrix*/
sequence_t *sdmodel_generate_sequences (sdmodel * mo, int seed,
                                        int global_len, long seq_number,
                                        int Tmax)
{
# define CUR_PROC "sdmodel_generate_sequences"

  /* An end state is characterized by not having out-going transition. */
  unsigned long tm;             /* Time seed */
  sequence_t *sq = NULL;
  int state, n, i, j, m, reject_os, reject_tmax, badseq, trans_class;
  double p, sum;
  int len = global_len, up = 0, stillbadseq = 0, reject_os_tmp = 0;
  int obsLength = 0;
  int silent_len = 0, badSilentStates = 0;
  int lastStateSilent = 0;
  int matchcount = 0;

  sq = sequence_calloc (seq_number);

  if (!sq) {
    mes_proc ();
    goto STOP;
  }
  if (len <= 0)
    /* A specific length of the sequences isn't given. As a model should have
       an end state, the konstant MAX_SEQ_LEN is used. */
    len = (int) MAX_SEQ_LEN;

  if (seed > 0) {
    GHMM_RNG_SET (RNG, seed);
  }
  else {
    tm = rand ();
    GHMM_RNG_SET (RNG, tm);
    fprintf (stderr, "# using rng '%s' seed=%ld\n", GHMM_RNG_NAME (RNG), tm);
  }

  n = 0;
  reject_os = reject_tmax = 0;

  while (n < seq_number) {
    /* Test: A new seed for each sequence */
    /*   ghmm_rng_timeseed(RNG); */
    stillbadseq = badseq = 0;
    matchcount = 0;
    ARRAY_CALLOC (sq->seq[n], len);
    ARRAY_CALLOC (sq->states[n], len);

    /* Get a random initial state i */
    p = GHMM_RNG_UNIFORM (RNG);
    sum = 0.0;
    for (i = 0; i < mo->N; i++) {
      sum += mo->s[i].pi;
      if (sum >= p)
        break;
    }
    /* assert( !mo->silent[i] ); */

    if (mo->model_type == kSilentStates) {

      if (!mo->silent[i]) {     /* fDte emits */
        lastStateSilent = 0;
        silent_len = 0;
        /* Get a random initial output m */
        p = GHMM_RNG_UNIFORM (RNG);
        sum = 0.0;
        for (m = 0; m < mo->M; m++) {
          sum += mo->s[i].b[m];
          if (sum >= p)
            break;
        }
        sq->seq[n][0] = m;
        sq->states[n][0] = i;
        if (mo->s[i].countme) {
          matchcount++;
        }
        state = 1;
      }
      else {                    /* silent state: we do nothing, no output */
        state = 0;
        lastStateSilent = 1;
        silent_len = 1;
      }
    }
    else {
      /* Get a random initial output m */
      p = GHMM_RNG_UNIFORM (RNG);
      sum = 0.0;
      for (m = 0; m < mo->M; m++) {
        sum += mo->s[i].b[m];
        if (sum >= p)
          break;
      }
      sq->seq[n][0] = m;
      sq->states[n][0] = i;
      if (mo->s[i].countme) {
        matchcount++;
      }
      state = 1;
    }

    obsLength = 0;

    /* class NOT dependent on observable but on path!!!! */
    while (state < len && obsLength < Tmax) {
      /* Get a new state */

      if (mo->cos > 1)
        trans_class = mo->get_class (sq->seq[n], state);
      else
        trans_class = 0;
      p = GHMM_RNG_UNIFORM (RNG);
      sum = 0.0;
      for (j = 0; j < mo->s[i].out_states; j++) {
        sum += mo->s[i].out_a[trans_class][j];
        if (sum >= p)
          break;
      }

      if (sum == 0.0) {
        if (mo->s[i].out_states > 0) {
          /* Repudiate the sequence, if all smo->s[i].out_a[class][.] == 0,
             that is, class "class" isn't used in the original data:
             go out of the while-loop, n should not be counted. */
          /* printf("Zustand %d, class %d, len %d out_states %d \n", i, class,
             state, smo->s[i].out_states); */
          badseq = 1;
          /* break; */

          /* Try: If the class is "empty", try the neighbour class;
             first, sweep down to zero; if still no success, sweep up to
             COS - 1. If still no success --> Repudiate the sequence. */
          if (trans_class > 0 && up == 0) {
            /*trans_class--;*/
            continue;
          }
          else {
            if (trans_class < mo->cos - 1) {
              /*trans_class++;          // as it is not dependent on time!*/
              up = 1;
              continue;
            }
            else {
              stillbadseq = 1;
              break;
            }
          }
        }
        else {
          /* An end state is reached, get out of the while-loop */
          break;
        }
      }

      i = mo->s[i].out_id[j];

      if (mo->model_type == kSilentStates && mo->silent[i]) {   /* Get a silent state i */
        silent_len++;
        if (silent_len >= Tmax) {
          badSilentStates = 1;
          break;                /* reject this sequence */
        }
        else {
          badSilentStates = 0;
          lastStateSilent = 1;
        }
      }
      else {
        lastStateSilent = 0;
        silent_len = 0;
        /* Get a random output m from state i */
        p = GHMM_RNG_UNIFORM (RNG);
        sum = 0.0;
        for (m = 0; m < mo->M; m++) {
          sum += mo->s[i].b[m];
          if (sum >= p)
            break;
        }
        sq->seq[n][state] = m;
        sq->states[n][state] = i;
        if (mo->s[i].countme) {
          matchcount++;
        }
        state++;
      }

      /* Decide the class for the next step */
      /*trans_class = mo->get_class(sq->seq[n],state);*/
      up = 0;
      obsLength++;
    }                           /* while (state < len) , global_len depends on the data */

    if (badseq) {
      reject_os_tmp++;
    }

    if (stillbadseq) {
      reject_os++;
      m_free (sq->seq[n]);
      m_free (sq->states[n]);
      /*  printf("cl %d, s %d, %d\n", class, i, n); */
    }
    else {
      if (badSilentStates) {
        reject_tmax++;
        m_free (sq->seq[n]);
        m_free (sq->states[n]);
      }
      else {
        if (obsLength > Tmax) {
          reject_tmax++;
          m_free (sq->seq[n]);
          m_free (sq->seq[n]);
        }
        else {
          if (state < len) {
            ARRAY_REALLOC (sq->seq[n], state);
            ARRAY_REALLOC (sq->states[n], state);
          }
          sq->seq_len[n] = state;
          /* sq->seq_label[n] = label; */
          /* vector_d_print(stdout, sq->seq[n], sq->seq_len[n]," "," ",""); */
          n++;
        }
      }
    }


    /*    printf("reject_os %d, reject_tmax %d\n", reject_os, reject_tmax); */
    if (reject_os > 10000) {
      mes_prot ("Reached max. no. of rejections\n");
      break;
    }
    if (!(n % 1000))
      printf ("%d Seqs. generated\n", n);
  }                             /* n-loop, while (n < seq_number) */


  if (reject_os > 0)
    printf ("%d sequences rejected (os)!\n", reject_os);
  if (reject_os_tmp > 0)
    printf ("%d sequences changed class\n", reject_os_tmp - reject_os);
  if (reject_tmax > 0)
    printf ("%d sequences rejected (Tmax, %d)!\n", reject_tmax, Tmax);


  return (sq);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */

  sequence_free (&sq);
  return (NULL);
# undef CUR_PROC
}                               /* data */

/*=======================================================================
generate sequences by calling generate_sequences_ext
========================================*/

/*
sequence_t *sdmodel_generate_sequences(sdmodel* mo, int seed, int global_len,
				     long seq_number, int Tmax) {
  sequence_ext_t *sq = sdmodel_generate_sequences_ext(mo,seed,global_len,seq_number,Tmax);
  printf("jaja\n\n");
  return sequence_unext(&sq);
}
*/
/*============================================================================*/
/* Some outputs */
/*============================================================================*/

void sdmodel_states_print (FILE * file, sdmodel * mo)
{
  int i, j;
  fprintf (file, "Modelparameters: \n M = %d \t N = %d\n", mo->M, mo->N);
  for (i = 0; i < mo->N; i++) {
    fprintf (file,
             "\nState %d \n PI = %.3f \n out_states = %d \n in_states = %d \n",
             i, mo->s[i].pi, mo->s[i].out_states, mo->s[i].in_states);
    fprintf (file, " Output probability:\t");
    for (j = 0; j < mo->M; j++)
      fprintf (file, "%.3f \t", mo->s[i].b[j]);
    fprintf (file, "\n Transition probability \n");
    fprintf (file, "  Out states (Id, a):\t");
    for (j = 0; j < mo->s[i].out_states; j++)
      fprintf (file, "FIXME: out_a is a matrix"/*(%d, %.3f) \t", mo->s[i].out_id[j], mo->s[i].out_a[j]*/);
    fprintf (file, "\n");
    fprintf (file, "  In states (Id, a):\t");
    for (j = 0; j < mo->s[i].in_states; j++)
      fprintf (file, "FIXME: in_a is a matrix"/*(%d, %.3f) \t", mo->s[i].in_id[j], mo->s[i].in_a[j]*/);
    fprintf (file, "\n");
  }
}                               /* model_states_print */

/*============================================================================*/

void sdmodel_Ak_print (FILE * file, sdmodel * mo, int k, char *tab,
                       char *separator, char *ending)
{
  int i, j, out_state;
  for (i = 0; i < mo->N; i++) {
    out_state = 0;
    fprintf (file, "%s", tab);
    if (mo->s[i].out_states > 0 && mo->s[i].out_id[out_state] == 0) {
      fprintf (file, "%.2f", mo->s[i].out_a[k][out_state]);
      out_state++;
    }
    else
      fprintf (file, "0.00");
    for (j = 1; j < mo->N; j++) {
      if (mo->s[i].out_states > out_state && mo->s[i].out_id[out_state] == j) {
        fprintf (file, "%s %.2f", separator, mo->s[i].out_a[k][out_state]);
        out_state++;
      }
      else
        fprintf (file, "%s 0.00", separator);
    }
    fprintf (file, "%s\n", ending);
  }
}                               /* model_A_print */

/*============================================================================*/

void sdmodel_B_print (FILE * file, sdmodel * mo, char *tab, char *separator,
                      char *ending)
{
  int i, j;
  for (i = 0; i < mo->N; i++) {
    fprintf (file, "%s", tab);
    fprintf (file, "%.2f", mo->s[i].b[0]);
    for (j = 1; j < mo->M; j++)
      fprintf (file, "%s %.2f", separator, mo->s[i].b[j]);
    fprintf (file, "%s\n", ending);
  }
}                               /* model_B_print */

/*============================================================================*/

void sdmodel_Pi_print (FILE * file, sdmodel * mo, char *tab, char *separator,
                       char *ending)
{
  int i;
  fprintf (file, "%s%.2f", tab, mo->s[0].pi);
  for (i = 1; i < mo->N; i++)
    fprintf (file, "%s %.2f", separator, mo->s[i].pi);
  fprintf (file, "%s\n", ending);
}                               /* model_Pi_print */


void model_to_sdmodel (const model * mo, sdmodel * smo, int klass)
{
#define CUR_PROC "model_to_sdmodel"
  int i, j, m, nachf, vorg;

  for (i = 0; i < mo->N; i++) {
    nachf = mo->s[i].out_states;
    vorg = mo->s[i].in_states;
    /* Copy the values */
    for (j = 0; j < nachf; j++) {
      smo->s[i].out_a[klass][j] = mo->s[i].out_a[j];
      smo->s[i].out_id[j] = mo->s[i].out_id[j];
    }
    for (j = 0; j < vorg; j++) {
      smo->s[i].in_a[klass][j] = mo->s[i].in_a[j];
      smo->s[i].in_id[j] = mo->s[i].in_id[j];
    }
    for (m = 0; m < mo->M; m++)
      smo->s[i].b[m] = mo->s[i].b[m];
    smo->s[i].pi = mo->s[i].pi;
    smo->s[i].out_states = nachf;
    smo->s[i].in_states = vorg;
  }
  smo->prior = mo->prior;
#undef CUR_PROC
}


model *sdmodel_to_model (const sdmodel * mo, int kclass)
{
#define CUR_PROC "sdmodel_to_model"
  /*
   * Set the pointer appropriately
   */
  int i, j, nachf, vorg, m;
  model *m2 = NULL;
  ARRAY_CALLOC (m2, 1);
  ARRAY_CALLOC (m2->s, mo->N);
  for (i = 0; i < mo->N; i++) {
    nachf = mo->s[i].out_states;
    vorg = mo->s[i].in_states;
    ARRAY_CALLOC (m2->s[i].out_id, nachf);
    ARRAY_CALLOC (m2->s[i].out_a, nachf);
    ARRAY_CALLOC (m2->s[i].in_id, vorg);
    ARRAY_CALLOC (m2->s[i].in_a, vorg);
    ARRAY_CALLOC (m2->s[i].b, mo->M);
    /* Copy the values */
    for (j = 0; j < nachf; j++) {
      m2->s[i].out_a[j] = mo->s[i].out_a[kclass][j];
      m2->s[i].out_id[j] = mo->s[i].out_id[j];
    }
    for (j = 0; j < vorg; j++) {
      m2->s[i].in_a[j] = mo->s[i].in_a[kclass][j];
      m2->s[i].in_id[j] = mo->s[i].in_id[j];
    }
    for (m = 0; m < mo->M; m++)
      m2->s[i].b[m] = mo->s[i].b[m];
    m2->s[i].pi = mo->s[i].pi;
    m2->s[i].out_states = nachf;
    m2->s[i].in_states = vorg;
  }
  m2->N = mo->N;
  m2->M = mo->M;
  m2->prior = mo->prior;

  m2->model_type = mo->model_type;
  if (mo->model_type == kSilentStates) {
    assert (mo->silent != NULL);
    ARRAY_CALLOC (m2->silent, mo->N);
    for (i = 0; i < m2->N; i++) {
      m2->silent[i] = mo->silent[i];
    }
    m2->topo_order_length = mo->topo_order_length;
    ARRAY_CALLOC (m2->topo_order, mo->topo_order_length);
    for (i = 0; i < m2->topo_order_length; i++) {
      m2->topo_order[i] = mo->topo_order[i];
    }
  }
  return (m2);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  m_free (m2->silent);
  m_free (m2->topo_order);
  model_free (&m2);
  return (NULL);
#undef CUR_PROC
}

/*============================================================================*/

/*============================================================================*/

/*state* state_copy(state *my_state) {
  state* new_state = (state*) malloc(sizeof(state));

  state_copy_to(my_state,new_state);

  return new_state;

  }*/ /* state_copy */

/*============================================================================*/

/*void state_copy_to(state *source, state* dest) {
  dest->pi         = source->pi;
  dest->out_states = source->out_states;
  dest->in_states  = source->in_states;
  dest->fix        = source->fix;

  dest->b          = malloc(xxx);
  memcpy(dest->b,source->b,xxx);

  dest->out_id     = malloc(xxx);
  memcpy(dest->out_id,source->out_id,xxx);

  dest->in_id      = malloc(xxx);
  memcpy(dest->in_id,source->in_id,xxx);

  dest->out_a      = malloc(xxx);
  memcpy(dest->out_a,source->out_a,xxx);

  dest->in_a       = malloc(xxx);
  memcpy(dest->in_a,source->in_a,xxx);
  }*/ /* state_copy_to */


/*===================== E n d   o f  f i l e  "model.c"       ===============*/
