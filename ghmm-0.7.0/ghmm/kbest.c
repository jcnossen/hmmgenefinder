/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/kbest.c
*       Authors:  Anyess von Bock, Alexander Riemer, Janne Grunau
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
*       This file is version $Revision: 1295 $
*                       from $Date: 2005-09-02 19:33:44 +0200 (Fri, 02 Sep 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "mes.h"
#include "mprintf.h"
#include "model.h"
#include "kbestbasics.h"
#include "kbest.h"
#include <ghmm/internal.h>


/**
  Builds logarithmic transition matrix from the states' in_a values
  the row for each state is the logarithmic version of the state's in_a
  @return transition matrix with logarithmic values, 1.0 if a[i,j] = 0
  @param s:           array of all states of the model
  @param N:           number of states in the model
 */
static double **kbest_buildLogMatrix (state * s, int N)
{
#define CUR_PROC "kbest_buildLogMatrix"
  int i, j;
  double **log_a;               /* log(a(i,j)) => log_a[i*N+j] */

  /* create & initialize matrix: */
  ARRAY_MALLOC (log_a, N);
  for (i = 0; i < N; i++) {
    ARRAY_MALLOC (log_a[i], s[i].in_states);
    for (j = 0; j < s[i].in_states; j++)
      log_a[i][j] = log (s[i].in_a[j]);
  }
  return log_a;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  mes_prot ("kbest_buildLogMatrix failed\n");
  exit (1);
#undef CUR_PROC
}


/**
   Calculates the most probable labeling for the given sequence in the given
   model using k-best decoding.
   Labels must be from interval [0:max_label] without gaps!!! (not checked)
   Model must not have silent states. (checked in Python wrapper)
   @return array of labels (internal representation)
   @param mo:         pointer to a model
   @param o_seq:      output sequence (array of internal representation chars)
   @param seq_len:    length of output sequence
   @param k:          number of hypotheses to keep for each state
   @param log_p:      variable reference to store the log prob. of the labeling
 */
int *kbest (model * mo, int *o_seq, int seq_len, int k, double *log_p)
{
#define CUR_PROC "kbest"
  int i, t, c, l, m;            /* counters */
  int no_oldHyps;               /* number of hypotheses until position t-1 */
  int b_index, i_id;            /* index for addressing states' b arrays */
  int no_labels = 0;
  int exists, g_nr;
  int *states_wlabel;
  int *label_max_out;
  char *str;

  /* logarithmized transition matrix A, log(a(i,j)) => log_a[i*N+j],
     1.0 for zero probability */
  double **log_a;

  /* matrix of hypotheses, holds for every position in the sequence a list
     of hypotheses */
  hypoList **h;
  hypoList *hP;

  /* vectors for rows in the matrices */
  int *hypothesis;

  /* pointer & prob. of the k most probable hypotheses for each state
     - matrices of dimensions #states x k:  argm(i,l) => argmaxs[i*k+l] */
  double *maxima;
  hypoList **argmaxs;

  /* pointer to & probability of most probable hypothesis in a certain state */
  hypoList *argmax;
  double sum;

  /* break if sequence empty or k<1: */
  if (seq_len <= 0 || k <= 0)
    return NULL;

  ARRAY_CALLOC (h, seq_len);

  /** 1. Initialization (extend empty hypothesis to #labels hypotheses of
         length 1): */

  /* get number of labels (= maximum label + 1)
     and number of states with those labels */
  ARRAY_CALLOC (states_wlabel, mo->N);
  ARRAY_CALLOC (label_max_out, mo->N);
  for (i = 0; i < mo->N; i++) {
    c = mo->s[i].label;
    states_wlabel[c]++;
    if (c > no_labels)
      no_labels = c;
    if (mo->s[i].out_states > label_max_out[c])
      label_max_out[c] = mo->s[i].out_states;
  }
  /* add one to the maximum label to get the number of labels */
  no_labels++;
  ARRAY_REALLOC (states_wlabel, no_labels);
  ARRAY_REALLOC (label_max_out, no_labels);

  /* initialize h: */
  hP = h[0];
  for (i = 0; i < mo->N; i++) {
    if (mo->s[i].pi > KBEST_EPS) {
      /* printf("Found State %d with initial probability %f\n", i, mo->s[i].pi); */
      exists = 0;
      while (hP != NULL) {
        if (hP->hyp_c == mo->s[i].label) {
          /* add entry to the gamma list */
          g_nr = hP->gamma_states;
          hP->gamma_id[g_nr] = i;
          hP->gamma_a[g_nr] =
            log (mo->s[i].pi) +
            log (mo->s[i].b[get_emission_index (mo, i, o_seq[0], 0)]);
          hP->gamma_states = g_nr + 1;
          exists = 1;
          break;
        }
        else
          hP = hP->next;
      }
      if (!exists) {
        hlist_insertElem (&(h[0]), mo->s[i].label, NULL);
        /* initiallize gamma-array with safe size (number of states) and add the first entry */
        ARRAY_MALLOC (h[0]->gamma_a, states_wlabel[mo->s[i].label]);
        ARRAY_MALLOC (h[0]->gamma_id, states_wlabel[mo->s[i].label]);
        h[0]->gamma_id[0] = i;
        h[0]->gamma_a[0] =
          log (mo->s[i].pi) +
          log (mo->s[i].b[get_emission_index (mo, i, o_seq[0], 0)]);
        h[0]->gamma_states = 1;
        h[0]->chosen = 1;
      }
      hP = h[0];
    }
  }
  /* reallocating the gamma list to the real size */
  hP = h[0];
  while (hP != NULL) {
    ARRAY_REALLOC (hP->gamma_a, hP->gamma_states);
    ARRAY_REALLOC (hP->gamma_id, hP->gamma_states);
    hP = hP->next;
  }

  /* calculate transition matrix with logarithmic values: */
  log_a = kbest_buildLogMatrix (mo->s, mo->N);

  /* initialize temporary arrays: */
  ARRAY_MALLOC (maxima, mo->N * k);                             /* for each state save k */
  ARRAY_MALLOC (argmaxs, mo->N * k);


  /*------ Main loop: Cycle through the sequence: ------*/
  for (t = 1; t < seq_len; t++) {

    /* put o_seq[t-1] in emission history: */
    update_emission_history (mo, o_seq[t - 1]);

    /** 2. Propagate hypotheses forward and update gamma: */
    no_oldHyps =
      hlist_propFwd (mo, h[t - 1], &(h[t]), no_labels, states_wlabel,
                     label_max_out);
    /* printf("t = %d (%d), no of old hypotheses = %d\n", t, seq_len, no_oldHyps); */

    /*-- calculate new gamma: --*/
    hP = h[t];
    /* cycle through list of hypotheses */
    while (hP != NULL) {

      for (i = 0; i < hP->gamma_states; i++) {
        /* if hypothesis hP ends with label of state i:
           gamma(i,c):= log(sum(exp(a(j,i)*exp(oldgamma(j,old_c)))))
           + log(b[i](o_seq[t]))
           else: gamma(i,c):= -INF (represented by 1.0) */
        i_id = hP->gamma_id[i];
        hP->gamma_a[i] = logGammaSum (log_a[i_id], &mo->s[i_id], hP->parent);
        b_index = get_emission_index (mo, i_id, o_seq[t], t);
        if (b_index < 0) {
          hP->gamma_a[i] = 1.0;
          if (mo->s[i_id].order > t)
            continue;
          else {
            str = mprintf (NULL, 0,
                           "i_id: %d, o_seq[%d]=%d\ninvalid emission index!\n",
                           i_id, t, o_seq[t]);
            mes_prot (str);
            m_free (str);
          }
        }
        else
          hP->gamma_a[i] += log (mo->s[i_id].b[b_index]);
        /*printf("%g = %g\n", log(mo->s[i_id].b[b_index]), hP->gamma_a[i]); */
        if (hP->gamma_a[i] > 0.0) {
          mes_prot ("gamma to large. kbest failed\n");
          exit (1);
        }
      }
      hP = hP->next;
    }

    /** 3. Choose the k most probable hypotheses for each state and discard all
	   hypotheses that were not chosen: */

    /* initialize temporary arrays: */
    for (i = 0; i < mo->N * k; i++) {
      maxima[i] = 1.0;
      argmaxs[i] = NULL;
    }

    /* cycle through hypotheses & calculate the k most probable hypotheses for
       each state: */
    hP = h[t];
    while (hP != NULL) {
      for (i = 0; i < hP->gamma_states; i++) {
        i_id = hP->gamma_id[i];
        if (hP->gamma_a[i] > KBEST_EPS)
          continue;
        /* find first best hypothesis that is worse than current hypothesis: */
        for (l = 0;
             l < k && maxima[i_id * k + l] < KBEST_EPS
             && maxima[i_id * k + l] > hP->gamma_a[i]; l++);
        if (l < k) {
          /* for each m>l: m'th best hypothesis becomes (m+1)'th best */
          for (m = k - 1; m > l; m--) {
            argmaxs[i_id * k + m] = argmaxs[i_id * k + m - 1];
            maxima[i_id * k + m] = maxima[i_id * k + m - 1];
          }
          /* save new l'th best hypothesis: */
          maxima[i_id * k + l] = hP->gamma_a[i];
          argmaxs[i_id * k + l] = hP;
        }
      }
      hP = hP->next;
    }

    /* set 'chosen' for all hypotheses from argmaxs array: */
    for (i = 0; i < mo->N * k; i++)
      /* only choose hypotheses whose prob. is at least threshold*max_prob */
      if (maxima[i] != 1.0
          && maxima[i] >= KBEST_THRESHOLD + maxima[(i % mo->N) * k])
        argmaxs[i]->chosen = 1;

    /* remove hypotheses that were not chosen from the lists: */
    /* remove all hypotheses till the first chosen one */
    while (h[t] != NULL && 0 == h[t]->chosen)
      hlist_removeElem (&(h[t]));
    /* remove all other not chosen hypotheses */
    if (!h[t]) {
      mes_prot ("No chosen hypothesis. kbest failed\n");
      exit (1);
    }
    hP = h[t];
    while (hP->next != NULL) {
      if (1 == hP->next->chosen)
        hP = hP->next;
      else
        hlist_removeElem (&(hP->next));
    }
  }
  /* dispose of temporary arrays: */
  free (states_wlabel);
  free (label_max_out);
  free (argmaxs);
  free (maxima);
  /* transition matrix is no longer needed from here on */
  free (log_a);

  /** 4. Save the hypothesis with the highest probability over all states: */
  hP = h[seq_len - 1];
  argmax = NULL;
  *log_p = 1.0;                 /* log_p will store log of maximum summed probability */
  while (hP != NULL) {
    /* sum probabilities for each hypothesis over all states: */
    sum = logSum (hP->gamma_a, hP->gamma_states);
    /* and select maximum sum */
    if (sum < KBEST_EPS && (*log_p == 1.0 || sum > *log_p)) {
      *log_p = sum;
      argmax = hP;
    }
    hP = hP->next;
  }

  /* found a valid path? */
  if (*log_p < KBEST_EPS) {
    /* yes: extract chosen hypothesis: */
    ARRAY_MALLOC (hypothesis, seq_len);
    for (i = seq_len - 1; i >= 0; i--) {
      hypothesis[i] = argmax->hyp_c;
      argmax = argmax->parent;
    }
  }
  else
    /* no: return 1.0 representing -INF and an empty hypothesis */
    hypothesis = NULL;

  /* dispose of calculation matrices: */
  hP = h[seq_len - 1];
  while (hP != NULL)
    hlist_removeElem (&hP);
  free (h);
  return hypothesis;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  mes_prot ("kbest failed\n");
  exit (1);
#undef CUR_PROC
}
