/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/sviterbi.c
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

/* #if __EXPERIMENTAL__ == 3 */

#include <float.h>
#include <math.h>
#include <ghmm/internal.h>
#include "sviterbi.h"
#include "matrix.h"
#include "smodel.h"
#include "mes.h"

typedef struct local_store_t {
  double **log_b;
  double *phi;
  double *phi_new;
  int **psi;
} local_store_t;

static local_store_t *sviterbi_alloc (smodel * smo, int T);
static int sviterbi_free (local_store_t ** v, int n, int T);

/*----------------------------------------------------------------------------*/
static local_store_t *sviterbi_alloc (smodel * smo, int T)
{
#define CUR_PROC "sviterbi_alloc"
  local_store_t *v = NULL;
  ARRAY_CALLOC (v, 1);
  v->log_b = stat_matrix_d_alloc (smo->N, T);
  if (!(v->log_b)) {
    mes_proc ();
    goto STOP;
  }
  ARRAY_CALLOC (v->phi, smo->N);
  ARRAY_CALLOC (v->phi_new, smo->N);
  v->psi = matrix_i_alloc (T, smo->N);
  if (!(v->psi)) {
    mes_proc ();
    goto STOP;
  }
  return (v);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sviterbi_free (&v, smo->N, T);
  return (NULL);
#undef CUR_PROC
}                               /* sviterbi_alloc */


/*----------------------------------------------------------------------------*/
static int sviterbi_free (local_store_t ** v, int n, int T)
{
#define CUR_PROC "sviterbi_free"
  mes_check_ptr (v, return (-1));
  if (!*v)
    return (0);
  stat_matrix_d_free (&((*v)->log_b));
  m_free ((*v)->phi);
  m_free ((*v)->phi_new);
  matrix_i_free (&((*v)->psi), T);
  m_free (*v);
  return (0);
#undef CUR_PROC
}                               /* sviterbi_free */


/*----------------------------------------------------------------------------*/

static void sviterbi_precompute (smodel * smo, double *O, int T,
                                 local_store_t * v)
{
#define CUR_PROC "sviterbi_precompute"
  int j, t, m;
  double cb;

  /* Precomputing of log(b_j(O_t)) */
  for (j = 0; j < smo->N; j++) {
    for (t = 0; t < T; t++) {
      cb = 0.0;
      for (m = 0; m < smo->M; m++)
        cb += smodel_calc_cmbm (smo, j, m, O[t]);
      if (cb == 0.0)            /* DBL_EPSILON ? */
        v->log_b[j][t] = -DBL_MAX;
      else
        v->log_b[j][t] = log (cb);
    }
  }
#undef CUR_PROC
}                               /* sviterbi_precompute */



/*============================================================================*/
int *sviterbi (smodel * smo, double *O, int T, double *log_p)
{
#define CUR_PROC "sviterbi"
  int *state_seq = NULL;
  int t, j, i, osc;
  double value, max_value;
  local_store_t *v;

  v = sviterbi_alloc (smo, T);
  if (!v) {
    mes_proc ();
    goto STOP;
  }
  ARRAY_CALLOC (state_seq, T);
  /* Precomputing of log(bj(ot)) */
  sviterbi_precompute (smo, O, T, v);

  /* Initialization for  t = 0 */
  for (j = 0; j < smo->N; j++) {
    if (smo->s[j].pi == 0.0 || v->log_b[j][0] == -DBL_MAX)
      v->phi[j] = -DBL_MAX;
    else
      v->phi[j] = log (smo->s[j].pi) + v->log_b[j][0];
  }

  /* Recursion */
  for (t = 1; t < T; t++) {
    if (smo->cos == 1) {
      osc = 0;
    }
    else {
      if (!smo->class_change->get_class) {
        printf ("ERROR: get_class not initialized\n");
        goto STOP;
      }
      /*printf("1: cos = %d, k = %d, t = %d\n",smo->cos,smo->class_change->k,t);*/
      osc =  smo->class_change->get_class (smo, O, smo->class_change->k, t - 1);
      if (osc >= smo->cos){
        printf("ERROR: get_class returned index %d but model has only %d classes !\n",osc,smo->cos);
        goto STOP;
      }
    }


    for (j = 0; j < smo->N; j++) {
      /* find maximum */
      /* max_phi = phi[i] + log_in_a[j][i] ... */
      max_value = -DBL_MAX;
      v->psi[t][j] = -1;
      for (i = 0; i < smo->s[j].in_states; i++) {
        if (v->phi[smo->s[j].in_id[i]] > -DBL_MAX &&
            log (smo->s[j].in_a[osc][i]) > -DBL_MAX) {
          value = v->phi[smo->s[j].in_id[i]] + log (smo->s[j].in_a[osc][i]);
          if (value > max_value) {
            max_value = value;
            v->psi[t][j] = smo->s[j].in_id[i];
          }
        }
      }
      /* no maximum found (state is never reached or b(O[t]) = 0 */
      if (max_value == -DBL_MAX || v->log_b[j][t] == -DBL_MAX)
        v->phi_new[j] = -DBL_MAX;
      else
        v->phi_new[j] = max_value + v->log_b[j][t];
    }                           /* for (j = 0; j < smo->N; j++) */

    /* replace old phi with new phi */
    for (j = 0; j < smo->N; j++)
      v->phi[j] = v->phi_new[j];
  }

  /* Termination */
  max_value = -DBL_MAX;
  state_seq[T - 1] = -1;
  for (j = 0; j < smo->N; j++)
    if (v->phi[j] != -DBL_MAX && v->phi[j] > max_value) {
      max_value = v->phi[j];
      state_seq[T - 1] = j;
    }
  if (max_value == -DBL_MAX) {
    /* sequence can't be build from model, no backtracking possible */
    *log_p = -DBL_MAX;
    mes_proc ();
    goto STOP;
    /*
       for (t = T - 2; t >= 0; t--)
       state_seq[t] = -1;    
     */
  }
  else {
    *log_p = max_value;
    /* Backtracking */
    for (t = T - 2; t >= 0; t--)
      state_seq[t] = v->psi[t + 1][state_seq[t + 1]];
  }
  sviterbi_free (&v, smo->N, T);
  return (state_seq);

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  /* Free the memory space... */
  sviterbi_free (&v, smo->N, T);
  m_free (state_seq);
  return NULL;
#undef CUR_PROC
}                               /* sviterbi */


/*  #endif */ /* __EXPERIMENTAL__ == 3 */
