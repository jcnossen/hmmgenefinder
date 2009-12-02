/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/model.c
*       Authors:  Benhard Knab, Bernd Wichern, Benjamin Georgi, Alexander Schliep
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
*       This file is version $Revision: 1526 $
*                       from $Date: 2005-12-19 14:25:30 +0100 (Mon, 19 Dec 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/


#ifdef WIN32
#  include "win_config.h"
#endif

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <float.h>
#include <math.h>
#include "model.h"
#include "matrix.h"
#include "sequence.h"
#include "const.h"
#include "rng.h"
#include "foba.h"
#include "mes.h"
#include "mprintf.h"
#include "string.h"
#include "ghmm.h"
#include "modelutil.h"
#include "vector.h"
#include <ghmm/internal.h>

#define  __EPS 10e-6

/* Using floating point exceptions */
#ifdef DEBUG_FPE
#include <fenv.h>
static void __attribute__ ((constructor)) trapfpe (void)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}
#endif

typedef enum DFSFLAG { DONE, NOTVISITED, VISITED } DFSFLAG;


/*typedef struct local_store_t {
  DFSFLAG *colors;
  int    *topo_order;
  int    topo_order_length;
} local_store_t;

static local_store_t *topo_alloc(model *mo, int len);
static int topo_free(local_store_t **v, int n, int cos, int len); */


/*----------------------------------------------------------------------------*/
 int model_ipow (const model * mo, int x, unsigned int n)
{
#define CUR_PROC "model_ipow"
  int result = 1;

  if (mo->pow_lookup && (mo->M == x) && ((int)n <= mo->maxorder + 1))
    return mo->pow_lookup[n];
  else {
    while (n != 0) {
      if (n & 1)
        result *= x;
      x *= x;
      n >>= 1;
    }
    return result;
  }
#undef CUR_PROC
}

/*----------------------------------------------------------------------------*/
static int model_state_alloc (state * s, int M, int in_states,
                              int out_states)
{
# define CUR_PROC "model_state_alloc"
  int res = -1;
  ARRAY_CALLOC (s->b, M);
  if (out_states > 0) {
    ARRAY_CALLOC (s->out_id, out_states);
    ARRAY_CALLOC (s->out_a, out_states);
  }
  if (in_states > 0) {
    ARRAY_CALLOC (s->in_id, in_states);
    ARRAY_CALLOC (s->in_a, in_states);
  }
  res = 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (res);
# undef CUR_PROC
}                               /* model_state_alloc */

/*----------------------------------------------------------------------------*/

static int model_copy_vectors (model * mo, int index, double **a_matrix,
                               double **b_matrix, double *pi, int *fix)
{
#define CUR_PROC "model_copy_vectors"
  int i, cnt_out = 0, cnt_in = 0;
  mo->s[index].pi = pi[index];
  mo->s[index].fix = fix[index];
  for (i = 0; i < mo->M; i++)
    mo->s[index].b[i] = b_matrix[index][i];
  for (i = 0; i < mo->N; i++) {
    if (a_matrix[index][i]) {   /* Transitions to a following state possible */
      if (cnt_out >= mo->s[index].out_states) {
        mes_proc ();
        return (-1);
      }
      mo->s[index].out_id[cnt_out] = i;
      mo->s[index].out_a[cnt_out] = a_matrix[index][i];
      cnt_out++;
    }
    if (a_matrix[i][index]) {   /* Transitions to a previous state possible */
      if (cnt_in >= mo->s[index].in_states) {
        mes_proc ();
        return (-1);
      }
      mo->s[index].in_id[cnt_in] = i;
      mo->s[index].in_a[cnt_in] = a_matrix[i][index];
      cnt_in++;
    }
  }
  return (0);
#undef CUR_PROC
}                               /* model_copy_vectors */



/*============================================================================*/

/* Old prototype:

model **model_read(char *filename, int *mo_number, int **seq,
			 const int *seq_len, int seq_number) { */

model **model_read (char *filename, int *mo_number)
{
#define CUR_PROC "model_read"
  int j;
  long new_models = 0;
  scanner_t *s = NULL;
  model **mo = NULL;
  *mo_number = 0;
  s = scanner_alloc (filename);
  if (!s) {
    mes_proc ();
    goto STOP;
  }
  while (!s->err && !s->eof) {
    scanner_get_name (s);
    scanner_consume (s, '=');
    if (s->err)
      goto STOP;
    if (!strcmp (s->id, "HMM") || !strcmp (s->id, "hmm")) {
      (*mo_number)++;
      /* more mem */
      ARRAY_REALLOC (mo, *mo_number);
      mo[*mo_number - 1] = model_direct_read (s, (int *) &new_models);
      if (!mo[*mo_number - 1]) {
        mes_proc ();
        goto STOP;
      }
      /* Copies the model, that has already been read. */
      if (new_models > 1) {
        /* "-1" because mo_number++ has already been done. */
        ARRAY_REALLOC (mo, *mo_number - 1 + new_models);
        for (j = 1; j < new_models; j++) {
          mo[*mo_number] = model_copy (mo[*mo_number - 1]);
          if (!mo[*mo_number]) {
            mes_proc ();
            goto STOP;
          }
          (*mo_number)++;
        }
      }
    }
    else if (!strcmp (s->id, "HMM_SEQ")) {
      model **tmp_mo = NULL;
      tmp_mo = model_from_sequence_ascii (s, &new_models);
      ARRAY_REALLOC (mo, (*mo_number + new_models));
      for (j = 0; j < new_models; j++) {
        if (!tmp_mo[j]) {
          mes_proc ();
          goto STOP;
        }
        mo[*mo_number] = tmp_mo[j];
        (*mo_number)++;
      }
    }
    else {
      scanner_error (s, "unknown identifier");
      goto STOP;
    }
    scanner_consume (s, ';');
    if (s->err)
      goto STOP;
  }                             /* while(!s->err && !s->eof) */
  
  scanner_free(&s);
  return mo;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  scanner_free(&s);
  return NULL;
#undef CUR_PROC
}                               /* model_read */



/*============================================================================*/

model *model_direct_read (scanner_t * s, int *multip)
{
#define CUR_PROC "model_direct_read"
  int i, m_read, n_read, a_read, b_read, pi_read, len, prior_read, fix_read;
  int mt_read=0;
  model *mo = NULL;
  double **a_matrix = NULL, **b_matrix = NULL;
  double *pi_vector = NULL;
  int *fix_vector = NULL;
  m_read = n_read = a_read = b_read = pi_read = prior_read = fix_read = 0;
  *multip = 1;                  /* default */
  ARRAY_CALLOC (mo, 1);

  scanner_consume (s, '{');
  if (s->err)
    goto STOP;
  while (!s->err && !s->eof && s->c - '}') {
    scanner_get_name (s);
    if (strcmp (s->id, "M") && strcmp (s->id, "N") && strcmp (s->id, "Pi")
        && strcmp (s->id, "A") && strcmp (s->id, "B") &&
        strcmp (s->id, "multip") && strcmp (s->id, "prior") &&
        strcmp (s->id, "fix_state") && strcmp(s->id, "ModelType")) {
      scanner_error (s, "unknown identifier");
      goto STOP;
    }
    scanner_consume (s, '=');
    if (s->err)
      goto STOP;
    if (!strcmp (s->id, "multip")) {
      *multip = scanner_get_int (s);
      if (*multip < 1) {        /* Doesn't make any sense! */
        *multip = 1;
        mes_prot ("Multip < 1 ignored\n");
      }
    }
    else if (!strcmp (s->id, "M")) {    /*Number of output values */
      if (m_read) {
        scanner_error (s, "identifier M twice");
        goto STOP;
      }
      mo->M = scanner_get_int (s);
      m_read = 1;
    }
    else if (!strcmp (s->id, "N")) {    /*Number of states */
      if (n_read) {
        scanner_error (s, "identifier N twice");
        goto STOP;
      }
      mo->N = scanner_get_int (s);
      ARRAY_CALLOC (mo->s, mo->N);
      n_read = 1;
    }
    else if (!strcmp (s->id, "A")) {    /*Transition probability */
      if (!n_read) {
        scanner_error (s, "need N as a range for A");
        goto STOP;
      }
      if (a_read) {
        scanner_error (s, "identifier A twice");
        goto STOP;
      }
      ARRAY_CALLOC (a_matrix, mo->N);
      scanner_get_name (s);
      if (!strcmp (s->id, "matrix")) {
        if (matrix_d_read (s, a_matrix, mo->N, mo->N)) {
          scanner_error (s, "unable to read matrix A");
          goto STOP;
        }
        a_read = 1;
      }
      else {
        scanner_error (s, "unknown identifier");
        goto STOP;
      }
    }
    else if (!strcmp (s->id, "B")) {    /*Output probability */
      if ((!n_read) || (!m_read)) {
        scanner_error (s, "need M and N as a range for B");
        goto STOP;
      }
      if (b_read) {
        scanner_error (s, "identifier B twice");
        goto STOP;
      }
      ARRAY_CALLOC (b_matrix, mo->N);
      scanner_get_name (s);
      if (!strcmp (s->id, "matrix")) {
        if (matrix_d_read (s, b_matrix, mo->N, mo->M)) {
          scanner_error (s, "unable to read matrix B");
          goto STOP;
        }
        b_read = 1;
      }
      else {
        scanner_error (s, "unknown identifier");
        goto STOP;
      }
    }
    else if (!strcmp (s->id, "prior")) {        /*A prior model */
      if (prior_read) {
        scanner_error (s, "identifier prior twice");
        goto STOP;
      }
      mo->prior = scanner_get_edouble (s);
      if ((mo->prior < 0 || mo->prior > 1) && mo->prior != -1) {
        scanner_error (s, "invalid model prior");
        goto STOP;
      }
      prior_read = 1;
    }
    else if (!strcmp (s->id, "ModelType")) {    /* Model type*/
      if (mt_read) {
        scanner_error(s, "identifier ModelType twice");
        goto STOP;
      }
      mo->model_type = scanner_get_int(s);
      if (mo->model_type & (kSilentStates + kLabeledStates + kHigherOrderEmissions)) {
	scanner_error(s, "unsupported Model Type");
	goto STOP;
      }
      mt_read = 1;
    }
    else if (!strcmp (s->id, "Pi")) {   /*Initial state probabilty */
      if (!n_read) {
        scanner_error (s, "need N as a range for Pi");
        goto STOP;
      }
      if (pi_read) {
        scanner_error (s, "identifier Pi twice");
        goto STOP;
      }
      scanner_get_name (s);
      if (!strcmp (s->id, "vector")) {
        scanner_consume (s, '{');
        if (s->err)
          goto STOP;
        pi_vector = scanner_get_double_earray (s, &len);
        if (len != mo->N) {
          scanner_error (s, "wrong number of elements in PI");
          goto STOP;
        }
        scanner_consume (s, ';');
        if (s->err)
          goto STOP;
        scanner_consume (s, '}');
        if (s->err)
          goto STOP;
        pi_read = 1;
      }
      else {
        scanner_error (s, "unknown identifier");
        goto STOP;
      }
    }
    else if (!strcmp (s->id, "fix_state")) {
      if (!n_read) {
        scanner_error (s, "need N as a range for fix_state");
        goto STOP;
      }
      if (fix_read) {
        scanner_error (s, "identifier fix_state twice");
        goto STOP;
      }
      scanner_get_name (s);
      if (!strcmp (s->id, "vector")) {
        scanner_consume (s, '{');
        if (s->err)
          goto STOP;
        fix_vector = scanner_get_int_array (s, &len);
        if (len != mo->N) {
          scanner_error (s, "wrong number of elements in fix_state");
          goto STOP;
        }
        scanner_consume (s, ';');
        if (s->err)
          goto STOP;
        scanner_consume (s, '}');
        if (s->err)
          goto STOP;
        fix_read = 1;
      }
      else {
        scanner_error (s, "unknown identifier");
        goto STOP;
      }
    }
    else {
      scanner_error (s, "unknown identifier");
      goto STOP;
    }
    scanner_consume (s, ';');
    if (s->err)
      goto STOP;
  }                             /* while(..s->c-'}') */
  scanner_consume (s, '}');
  if (s->err)
    goto STOP;

  /* No prior read --> give it the value -1 */
  if (prior_read == 0)
    mo->prior = -1;
  /* Allocate memory for the model */
  for (i = 0; i < mo->N; i++) {
    mo->s[i].out_states = matrix_d_notzero_columns (a_matrix, i, mo->N);
    mo->s[i].in_states = matrix_d_notzero_rows (a_matrix, i, mo->N);
    if (model_state_alloc (mo->s + i, mo->M, mo->s[i].in_states,
                           mo->s[i].out_states)) {
      mes_proc ();
      goto STOP;
    }

    /* Assign the parameters to the model */
    if (!a_matrix) {
      fprintf (stderr, "no A matrix specified in file!\n");
      exit (1);
    }
    if (!b_matrix) {
      fprintf (stderr, "no B matrix specified in file!\n");
      exit (1);
    }
    if (!fix_vector) {
      fprintf (stderr, "no fix_state vector specified in file!\n");
      exit (1);
    }
    if (!pi_vector) {
      fprintf (stderr, "no Pi vector specified in file!\n");
      exit (1);
    }

    if (model_copy_vectors (mo, i, a_matrix, b_matrix, pi_vector, fix_vector)) {
      mes_proc ();
      goto STOP;
    }
  }
  ARRAY_CALLOC(mo->silent, mo->N);

  matrix_d_free (&a_matrix, mo->N);
  matrix_d_free (&b_matrix, mo->N);
  m_free (pi_vector);
  return (mo);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  matrix_d_free (&a_matrix, mo->N);
  matrix_d_free (&b_matrix, mo->N);
  m_free (pi_vector);
  model_free (&mo);
  return NULL;
#undef CUR_PROC
}                               /* model_direct_read */

/*============================================================================*/
/* Produces models from given sequences */
model **model_from_sequence (sequence_t * sq, long *mo_number)
{
#define CUR_PROC "model_from_sequence"
  long i;
  int max_symb;
  model **mo = NULL;
  ARRAY_CALLOC (mo, sq->seq_number);
  max_symb = sequence_max_symbol (sq);
  for (i = 0; i < sq->seq_number; i++)
    mo[i] = model_generate_from_sequence (sq->seq[i], sq->seq_len[i],
                                          max_symb + 1);
  *mo_number = sq->seq_number;
  return mo;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  for (i = 0; i < *mo_number; i++)
    model_free (&(mo[i]));
  return NULL;
#undef CUR_PROC
}                               /* model_from_sequence */

/*============================================================================*/
/* Produces models form given sequences */
model **model_from_sequence_ascii (scanner_t * s, long *mo_number)
{
#define CUR_PROC "model_from_sequence_ascii"
  long i;
  int max_symb;
  model **mo = NULL;
  sequence_t *sq = NULL;

  scanner_consume (s, '{');
  if (s->err)
    goto STOP;
  while (!s->err && !s->eof && s->c - '}') {
    scanner_get_name (s);
    scanner_consume (s, '=');
    if (s->err)
      goto STOP;
    /* Reads sequences on normal format */
    if (!strcmp (s->id, "SEQ")) {
      sq = sequence_read_alloc (s);
      if (!sq) {
        mes_proc ();
        goto STOP;
      }
    }
    else {
      scanner_error (s, "unknown identifier");
      goto STOP;
    }
    scanner_consume (s, ';');
    if (s->err)
      goto STOP;
  }                             /* while(..s->c-'}') */
  scanner_consume (s, '}');
  if (s->err)
    goto STOP;

  ARRAY_CALLOC (mo, sq->seq_number);
  /* The biggest symbol that occurs */
  max_symb = sequence_max_symbol (sq);
  for (i = 0; i < sq->seq_number; i++)
    mo[i] = model_generate_from_sequence (sq->seq[i], sq->seq_len[i],
                                          max_symb + 1);

  *mo_number = sq->seq_number;
  sequence_free (&sq);
  return mo;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_free (&sq);
  for (i = 0; i < *mo_number; i++)
    model_free (&(mo[i]));
  return NULL;
#undef CUR_PROC
}                               /* model_from_sequence_ascii */

/*============================================================================*/

int model_free (model ** mo)
{
#define CUR_PROC "model_free"
  int i,j;
  mes_check_ptr (mo, return (-1));

  for (i = 0; i < (*mo)->N; i++)
    state_clean (&(*mo)->s[i]);

  if ((*mo)->s){
    m_free ((*mo)->s);
  }
  if ((*mo)->silent){
    m_free ((*mo)->silent);
  }
  if ((*mo)->tied_to) {
    m_free ((*mo)->tied_to);
  }
  if ((*mo)->topo_order){
    m_free ((*mo)->topo_order);
  }
  if ((*mo)->pow_lookup){
    m_free ((*mo)->pow_lookup);
  }
  
  

  /* Optional attributes for storing representation information from the XML */  
  if ((*mo)->alphabet){
    for(i=0; i < (*mo)->S; i++) {
      for(j=0; j < (*mo)->alphabet_size[i]; j++) {
	m_free ((*mo)->alphabet[i][j]);	
      }
    }
    m_free ((*mo)->alphabet);
    m_free ((*mo)->alphabet_size);
  }
  if ((*mo)->position){
    m_free ((*mo)->position);
  }
  if ((*mo)->label_alphabet){
    for(i=0; i < (*mo)->label_size; i++) {
      m_free ((*mo)->label_alphabet[i]);
    }
    m_free ((*mo)->label_alphabet); 
  }

  
  
  
  m_free (*mo);
  return (0);
#undef CUR_PROC
}  /* model_free */


int model_free_background_distributions (background_distributions * bg)
{
#define CUR_PROC "model_free_background_distributions"

  if (bg->order){
    m_free (bg->order);
  }
  if (bg->b){
    matrix_d_free (&(bg->b), bg->n);
  }
  m_free (bg);

  return (0);
#undef CUR_PROC
}



/*============================================================================*/
model *model_copy (const model * mo)
{
# define CUR_PROC "model_copy"
  int i, j, nachf, vorg, m, size;
  model *m2 = NULL;
  ARRAY_CALLOC (m2, 1);
  ARRAY_CALLOC (m2->s, mo->N);
  ARRAY_CALLOC (m2->silent, mo->N);
  if (mo->model_type & kTiedEmissions) {
    ARRAY_CALLOC (m2->tied_to, mo->N);
  }
  else
    m2->tied_to = NULL;

  if (mo->model_type & kHasBackgroundDistributions) {
    ARRAY_CALLOC (m2->background_id, mo->N);
    m2->bp = mo->bp;
  }
  else
    m2->background_id = NULL;
  
  ARRAY_MALLOC (m2->pow_lookup, mo->maxorder + 2);
  
  for (i = 0; i < mo->N; i++) {
    nachf = mo->s[i].out_states;
    vorg = mo->s[i].in_states;
    
    ARRAY_CALLOC (m2->s[i].out_id, nachf);
    ARRAY_CALLOC (m2->s[i].out_a, nachf);
    ARRAY_CALLOC (m2->s[i].in_id, vorg);
    ARRAY_CALLOC (m2->s[i].in_a, vorg);
    /* allocate enough memory for higher order states */
    ARRAY_CALLOC (m2->s[i].b, model_ipow (mo, mo->M, mo->s[i].order + 1));

    /* copy the values */
    for (j = 0; j < nachf; j++) {
      m2->s[i].out_a[j] = mo->s[i].out_a[j];
      m2->s[i].out_id[j] = mo->s[i].out_id[j];
    }
    for (j = 0; j < vorg; j++) {
      m2->s[i].in_a[j] = mo->s[i].in_a[j];
      m2->s[i].in_id[j] = mo->s[i].in_id[j];
    }
    /* copy all b values for higher order states */
    size = model_ipow (mo, mo->M, mo->s[i].order + 1);
    for (m = 0; m < size; m++)
      m2->s[i].b[m] = mo->s[i].b[m];

    m2->s[i].pi = mo->s[i].pi;
    m2->s[i].fix = mo->s[i].fix;
    if (mo->model_type & kSilentStates)
      m2->silent[i] = mo->silent[i];
    if (mo->model_type & kTiedEmissions)
      m2->tied_to[i] = mo->tied_to[i];
    if (mo->model_type & kLabeledStates)
      m2->s[i].label = mo->s[i].label;
    if (mo->model_type & kHigherOrderEmissions)
      m2->s[i].order = mo->s[i].order;
    if (mo->model_type & kHasBackgroundDistributions)
      m2->background_id[i] = mo->background_id[i];
    m2->s[i].out_states = nachf;
    m2->s[i].in_states = vorg;
  }

  m2->N = mo->N;
  m2->M = mo->M;
  m2->prior = mo->prior;
  
  for (i = 0; i < mo->maxorder + 2; i++){
    m2->pow_lookup[i] = mo->pow_lookup[i];
  }
  
  m2->model_type = mo->model_type;
  /* not necessary but the history is at least initialised */
  m2->emission_history = mo->emission_history;
  return (m2);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  model_free (&m2);
  return (NULL);
# undef CUR_PROC
}                               /* model_copy */


/*============================================================================*/
int model_check (const model * mo)
{
# define CUR_PROC "model_check"
  int res = -1;
  double sum;
  int i, j, imag=0;
  char *str;
  /* The sum of the Pi[i]'s is 1 */

  sum = 0.0;
  for (i = 0; i < mo->N; i++) {
    sum += mo->s[i].pi;
  }

  if (fabs (sum - 1.0) >= EPS_PREC) {
    mes_prot ("sum Pi[i] != 1.0\n");
    goto STOP;
  }

  /* check each state */
  for (i = 0; i < mo->N; i++) {
    sum = 0.0;
    if (mo->s[i].out_states == 0) {
      str = mprintf (NULL, 0, "out_states = 0 (state %d -> final state!)\n", i);
      mes_prot (str);
    }
    /* Sum the a[i][j]'s : normalized out transitions */
    for (j = 0; j < mo->s[i].out_states; j++) {
      sum += mo->s[i].out_a[j];
      /* printf ("    out_a[%d][%d] = %8.5f\n", i, j, mo->s[i].out_a[j]); */
    }
    if (fabs (sum - 1.0) >= EPS_PREC) {
      str = mprintf (NULL, 0, "sum out_a[j] = %.5f != 1.0 (state %d)\n",
                           sum, i);
      mes_prot (str);
      m_free (str);
      /* goto STOP; */
    }
    /* Sum the a[i][j]'s : normalized in transitions */
    sum = mo->s[i].pi;
    for (j=0; j<mo->s[i].in_states; j++) {
      sum += mo->s[i].in_a[j];
    }
    if (fabs (sum) == 0.0) {
      imag = 1;
      str = mprintf (NULL, 0, "state %d can't be reached\n", i);
      mes_prot (str);
      m_free (str);
    }
    else if (fabs (sum-1.0) >= EPS_PREC) {
      str = mprintf (NULL, 0, "sum out_a[j] = %.5f != 1.0 (state %d)\n", sum, i);
      mes_prot (str);
      m_free (str);
    }
    /* Sum the b[j]'s: normalized emission probs */
    sum = 0.0;
    for (j = 0; j < mo->M; j++)
      sum += mo->s[i].b[j];
    /* silent states */
    if ((mo->model_type & kSilentStates) && mo->silent[i] && (sum != 0.0))
      goto STOP;
    /* not reachable states */
    else if (imag && (fabs (sum + mo->M) >= EPS_PREC))
      goto STOP;
    /* normal states */
    if (fabs (sum - 1.0) >= EPS_PREC) {
      str = mprintf (NULL, 0, "sum b[j] = %.2f != 1.0 (state %d)\n", sum, i);
      mes_prot (str);
      m_free (str);
      goto STOP;
    }                           /* i over all states */
  }

  res = 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (res);
# undef CUR_PROC
}                               /* model_check */

/*============================================================================*/
int model_check_compatibility (model ** mo, int model_number)
{
#define CUR_PROC "model_check_compatibility"
  int i;
  for (i = 1; i < model_number; i++)
    if (-1 == model_check_compatibel_models (mo[0], mo[i]))
      return -1;

  return 0;
#undef CUR_PROC
}

/*============================================================================*/
int model_check_compatibel_models (const model * mo, const model * m2)
{
#define CUR_PROC "model_check_compatibel_models"
  int i, j;
  char *str;

  if (mo->N != m2->N) {
    str = mprintf(NULL, 0, "ERROR: different number of states (%d != %d)\n",
                   mo->N, m2->N);
    goto STOP;
  }
  if (mo->M != m2->M) {
    str = mprintf(NULL, 0, "ERROR: different number of possible outputs (%d != %d)\n",
                   mo->M, m2->M);
    goto STOP;
  }
  for (i=0; i<mo->N; ++i) {
    if (mo->s[i].out_states != m2->s[i].out_states) {
      str = mprintf(NULL, 0, "ERROR: different number of outstates (%d != %d) in state %d.\n",
                   mo->s[i].out_states, m2->s[i].out_states, i);
      goto STOP;
    }
    for (j=0; j<mo->s[i].out_states; ++j) {
      if (mo->s[i].out_id[j] != m2->s[i].out_id[j]) {
	str = mprintf(NULL, 0, "ERROR: different out_ids (%d != %d) in entry %d of state %d.\n",
		      mo->s[i].out_id[j], m2->s[i].out_id[j], j, i);
	goto STOP;
      }
    }
  }

  return 0;
STOP:
  mes_prot (str);
  m_free (str);
  return (-1);
#undef CUR_PROC
}                               /* model_check_compatibility */

/*============================================================================*/
model *model_generate_from_sequence (const int *seq, int seq_len,
                                     int anz_symb)
{
#define CUR_PROC "model_generate_from_sequence"
  int i;
  model *mo = NULL;
  state *s = NULL;
  ARRAY_CALLOC (mo, 1);
  mo->N = seq_len;
  mo->M = anz_symb;
  /* All models generated from sequences have to be LeftRight-models */
  mo->model_type = kLeftRight;

  /* Allocate memory for all vectors */
  ARRAY_CALLOC (mo->s, mo->N);
  for (i = 0; i < mo->N; i++) {
    if (i == 0) {               /* Initial state */
      if (model_state_alloc (mo->s, mo->M, 0, 1)) {
        mes_proc ();
        goto STOP;
      }
    }
    else if (i == mo->N - 1) {  /* End state */
      if (model_state_alloc (mo->s + i, mo->M, 1, 0)) {
        mes_proc ();
        goto STOP;
      }
    }
    else {                      /* others */
      if (model_state_alloc (mo->s + i, mo->M, 1, 1)) {
        mes_proc ();
        goto STOP;
      }
    }
  }

  /* Allocate states with the right values, the initial state and the end 
     state extra */
  for (i = 1; i < mo->N - 1; i++) {
    s = mo->s + i;
    s->pi = 0.0;
    s->out_states = 1;
    s->in_states = 1;
    s->b[seq[i]] = 1.0;         /* others stay 0 */
    *(s->out_id) = i + 1;
    *(s->in_id) = i - 1;
    *(s->out_a) = *(s->in_a) = 1.0;
  }

  /* Initial state */
  s = mo->s;
  s->pi = 1.0;
  s->out_states = 1;
  s->in_states = 0;
  s->b[seq[0]] = 1.0;
  *(s->out_id) = 1;
  *(s->out_a) = 1.0;
  /* No in_id and in_a */

  /* End state */
  s = mo->s + mo->N - 1;
  s->pi = 0.0;
  s->out_states = 0;
  s->in_states = 1;
  s->b[seq[mo->N - 1]] = 1.0;   /* All other b's stay zero */
  *(s->in_id) = mo->N - 2;
  *(s->in_a) = 1.0;
  /* No out_id and out_a */

  if (model_check (mo)) {
    mes_proc ();
    goto STOP;
  }
  return mo;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  model_free (&mo);
  return NULL;
#undef CUR_PROC
}                               /* model_generate_from_sequence */


/*===========================================================================*/

 static int get_random_output (model * mo, int i, int position)
{
#define CUR_PROC "get_random_output"
  int m, e_index;
  double p, sum=0.0;

  p = GHMM_RNG_UNIFORM (RNG);

  for (m = 0; m < mo->M; m++) {
    /* get the right index for higher order emission models */
    e_index = get_emission_index(mo, i, m, position);

    /* get the probability, exit, if the index is -1 */
    if (-1 != e_index) {
      sum += mo->s[i].b[e_index];
      if (sum >= p)
        break;
    }
    else {
      fprintf (stderr,
               "ERROR: State has order %d, but in the history are only %d emissions.\n",
               mo->s[i].order, position);
      return -1;
    }
  }

  if (mo->M == m) {
    fprintf (stderr,
             "ERROR: no valid output choosen. Are the Probabilities correct? sum: %g, p: %g\n",
             sum, p);
    return -1;
  }

  return (m);
#undef CUR_PROC
}                               /* get_random_output */

/*============================================================================*/

sequence_t *model_generate_sequences (model * mo, int seed, int global_len,
                                      long seq_number, int Tmax)
{
# define CUR_PROC "model_generate_sequences"

  /* An end state is characterized by not having an output probabiliy. */

  sequence_t *sq = NULL;
  int i, j, m;
  double p, sum;
  int len = global_len;
  /* int silent_len = 0; */
  int n = 0;
  int state = 0;

  sq = sequence_calloc (seq_number);
  if (!sq) {mes_proc (); goto STOP;}
  if (len <= 0)
    /* A specific length of the sequences isn't given. As a model should have
       an end state, the konstant MAX_SEQ_LEN is used. */
    len = (int) MAX_SEQ_LEN;

  if (seed > 0) {
    GHMM_RNG_SET (RNG, seed);
  }

  /* initialize the emission history */
  mo->emission_history = 0;

  while (n < seq_number) {
    /* printf("sequenz n = %d\n",n); */

    ARRAY_CALLOC (sq->seq[n], len);
    state = 0;
    
    /* Get a random initial state i */
    p = GHMM_RNG_UNIFORM (RNG);
    sum = 0.0;
    for (i = 0; i < mo->N; i++) {
      sum += mo->s[i].pi;
      if (sum >= p)
        break;
    }

    if ((mo->model_type & kHigherOrderEmissions) && (mo->s[i].order > 0)) {
      fprintf (stderr,
               "ERROR: State %d has emission order %d, but it's initial probability is not 0.\n",
               i, mo->s[i].order);
      exit (1);
    }

    if (mo->model_type & kSilentStates && mo->silent[i]) {
      /* silent state: we do nothing, no output */
      /* printf("first state %d silent\n",i);
         state = 0; 
         silent_len = silent_len + 1; */
    }
    else {
      /* first state emits */
      /* printf("first state %d not silent\n",i); */
      /* Get a random initial output m */
      m = get_random_output (mo, i, state);
      update_emission_history (mo, m);
      sq->seq[n][state++] = m;
    }

    /* check whether sequence was completed by inital state */
    if (state >= len) {
      /* printf("assinging length %d to sequence %d\n",state,n);
         printf("sequence complete...\n"); */
      sq->seq_len[n++] = state;
      continue;
    }

    while (state < len) {
      /* Get a random state i */
      p = GHMM_RNG_UNIFORM (RNG);
      sum = 0.0;
      for (j = 0; j < mo->s[i].out_states; j++) {
        sum += mo->s[i].out_a[j];
        if (sum >= p)
          break;
      }

      /* printf("state %d selected (i: %d, j: %d) at position %d\n",mo->s[i].out_id[j],i,j,state); */

      if (sum == 0.0) {
	/* printf("final state reached - aborting\n"); */
	/* Set Sequence length and sample the next sequence */
	sq->seq_len[n++] = state;
	break;
      }

      i = mo->s[i].out_id[j];
      if ((mo->model_type & kSilentStates) && mo->silent[i]) {  /* Get a silent state i */
        /* printf("silent state \n");
           silent_len += 1; */
        /*if (silent_len >= Tmax) {
           printf("%d silent states reached -> silent circle - aborting...\n",silent_len);
           sq->seq_len[n++] = state; 
           break;
           } */
      }
      else {
        /* Get a random output m from state i */
        m = get_random_output (mo, i, state);
        update_emission_history (mo, m);
        sq->seq[n][state++] = m;
      }
      if (state == len)
        sq->seq_len[n++] = state;

    }                           /* while (state < len) */
  }                             /* while( n < seq_number ) */

  return (sq);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_free (&sq);
  return (NULL);
# undef CUR_PROC
}                               /* data */

/*============================================================================*/

double model_likelihood (model * mo, sequence_t * sq)
{
# define CUR_PROC "model_likelihood"
  double log_p_i, log_p;
  int found, i;
  char *str;

  /* printf("***  model_likelihood:\n"); */

  found = 0;
  log_p = 0.0;
  for (i = 0; i < sq->seq_number; i++) {

/* 	printf("sequence:\n"); */
/* 	for (j=0;j < sq->seq_len[i];j++) {  */
/* 		printf("%d, ",sq->seq[i][j]); */
/* 	} */
/* 	printf("\n"); */


    if (foba_logp (mo, sq->seq[i], sq->seq_len[i], &log_p_i) == -1) {
      mes_proc ();
      goto STOP;
    }

/* 	printf("\nlog_p_i = %f\n", log_p_i); */

    if (log_p_i != +1) {
      log_p += log_p_i;
      found = 1;
    }
    else {
      str = mprintf (NULL, 0, "sequence[%d] can't be build.\n", i);
      mes_prot (str);
    }
  }
  if (!found)
    log_p = +1.0;
  return (log_p);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return -1;
# undef CUR_PROC
}                               /* model_likelihood */



void model_set_transition (model * mo, int i, int j, double prob)
{
# define CUR_PROC "model_set_transition"
  int in, out;

  if (mo->s && mo->s[i].out_a && mo->s[j].in_a) {
    for (out = 0; out < mo->s[i].out_states; out++) {
      if (mo->s[i].out_id[out] == j) {
        mo->s[i].out_a[out] = prob;
        fprintf (stderr, "model_set_transition(0):State %d, %d, = %f\n", i, j,
                 prob);
        break;
      }
    }

    for (in = 0; in < mo->s[j].in_states; in++) {
      if (mo->s[j].in_id[in] == i) {
        mo->s[j].in_a[in] = prob;
        break;
      }
    }
  }
# undef CUR_PROC
}

/* model_set_transition */




/*============================================================================*/
/* Some outputs */
/*============================================================================*/

void model_states_print (FILE * file, model * mo)
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
      fprintf (file, "(%d, %.3f) \t", mo->s[i].out_id[j], mo->s[i].out_a[j]);
    fprintf (file, "\n");
    fprintf (file, "  In states (Id, a):\t");
    for (j = 0; j < mo->s[i].in_states; j++)
      fprintf (file, "(%d, %.3f) \t", mo->s[i].in_id[j], mo->s[i].in_a[j]);
    fprintf (file, "\n");
  }
}                               /* model_states_print */

/*============================================================================*/

void model_A_print (FILE * file, model * mo, char *tab, char *separator,
                    char *ending)
{
  int i, j, out_state;
  for (i = 0; i < mo->N; i++) {
    out_state = 0;
    fprintf (file, "%s", tab);
    if (mo->s[i].out_states > 0 && mo->s[i].out_id[out_state] == 0) {
      fprintf (file, "%.2f", mo->s[i].out_a[out_state]);
      out_state++;
    }
    else
      fprintf (file, "0.00");
    for (j = 1; j < mo->N; j++) {
      if (mo->s[i].out_states > out_state && mo->s[i].out_id[out_state] == j) {
        fprintf (file, "%s %.2f", separator, mo->s[i].out_a[out_state]);
        out_state++;
      }
      else
        fprintf (file, "%s 0.00", separator);
    }
    fprintf (file, "%s\n", ending);
  }
}                               /* model_A_print */

/*============================================================================*/

void model_B_print (FILE * file, model * mo, char *tab, char *separator,
                    char *ending)
{
  int i, j, size;

  for (i = 0; i < mo->N; i++) {
    fprintf (file, "%s", tab);
    fprintf (file, "%.2f", mo->s[i].b[0]);
    if (!(mo->model_type & kHigherOrderEmissions)) {
      for (j = 1; j < mo->M; j++)
        fprintf (file, "%s %.2f", separator, mo->s[i].b[j]);
      fprintf (file, "%s\n", ending);
    }
    else {
      size = model_ipow (mo, mo->M, mo->s[i].order + 1);
      for (j = 1; j < size; j++)
        fprintf (file, "%s %.2f", separator, mo->s[i].b[j]);
      fprintf (file, "%s\n", ending);
    }
  }
}                               /* model_B_print */

/*============================================================================*/

void model_Pi_print (FILE * file, model * mo, char *tab, char *separator,
                     char *ending)
{
  int i;
  fprintf (file, "%s%.2f", tab, mo->s[0].pi);
  for (i = 1; i < mo->N; i++)
    fprintf (file, "%s %.2f", separator, mo->s[i].pi);
  fprintf (file, "%s\n", ending);
}                               /* model_Pi_print */

void model_fix_print (FILE * file, model * mo, char *tab, char *separator,
                      char *ending)
{
  int i;
  fprintf (file, "%s%d", tab, mo->s[0].fix);
  for (i = 1; i < mo->N; i++)
    fprintf (file, "%s %d", separator, mo->s[i].fix);
  fprintf (file, "%s\n", ending);
}                               /* model_Pi_print */

/*============================================================================*/

void model_A_print_transp (FILE * file, model * mo, char *tab,
                           char *separator, char *ending)
{
# define CUR_PROC "model_A_print_transp"
  int i, j;
  int *out_state;

  ARRAY_CALLOC (out_state, mo->N);
  for (i = 0; i < mo->N; i++)
    out_state[i] = 0;

  for (j = 0; j < mo->N; j++) {
    fprintf (file, "%s", tab);
    if (mo->s[0].out_states != 0 && mo->s[0].out_id[out_state[0]] == j) {
      fprintf (file, "%.2f", mo->s[0].out_a[out_state[0]]);
      (out_state[0])++;
    }
    else
      fprintf (file, "0.00");
    for (i = 1; i < mo->N; i++) {
      if (mo->s[i].out_states != 0 && mo->s[i].out_id[out_state[i]] == j) {
        fprintf (file, "%s %.2f", separator, mo->s[i].out_a[out_state[i]]);
        (out_state[i])++;
      }
      else
        fprintf (file, "%s 0.00", separator);
    }
    fprintf (file, "%s\n", ending);
  }
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  m_free (out_state);
  return;
# undef CUR_PROC
}                               /* model_A_print_transp */

/*============================================================================*/

void model_B_print_transp (FILE * file, model * mo, char *tab,
                           char *separator, char *ending)
{
  int i, j;
  for (j = 0; j < mo->M; j++) {
    fprintf (file, "%s", tab);
    fprintf (file, "%.2f", mo->s[0].b[j]);
    for (i = 1; i < mo->N; i++)
      fprintf (file, "%s %.2f", separator, mo->s[i].b[j]);
    fprintf (file, "%s\n", ending);
  }
}                               /* model_B_print_transp */

/*============================================================================*/

void model_Pi_print_transp (FILE * file, model * mo, char *tab, char *ending)
{
  int i;
  for (i = 0; i < mo->N; i++)
    fprintf (file, "%s%.2f%s\n", tab, mo->s[i].pi, ending);
}                               /* model_Pi_print_transp */

/*============================================================================*/

void model_label_print (FILE * file, model * mo, char *tab, char *separator,
                        char *ending)
{
  int i;
  fprintf (file, "%s%d", tab, mo->s[0].label);
  for (i = 1; i < mo->N; i++)
    fprintf (file, "%s %d", separator, mo->s[i].label);
  fprintf (file, "%s\n", ending);
}                               /* model_label_print */

/*============================================================================*/
void model_print (FILE * file, model * mo)
{
  fprintf (file, "HMM = {\n\tM = %d;\n\tN = %d;\n", mo->M, mo->N);
  fprintf (file, "\tprior = %.3f;\n", mo->prior);
  fprintf (file, "\tModelType = %d;\n", mo->model_type);
  fprintf (file, "\tA = matrix {\n");
  model_A_print (file, mo, "\t", ",", ";");
  fprintf (file, "\t};\n\tB = matrix {\n");
  model_B_print (file, mo, "\t", ",", ";");
  fprintf (file, "\t};\n\tPi = vector {\n");
  model_Pi_print (file, mo, "\t", ",", ";");
  fprintf (file, "\t};\n\tfix_state = vector {\n");
  model_fix_print (file, mo, "\t", ",", ";");
  if (mo->model_type & kLabeledStates) {
    fprintf (file, "\t};\n\tlabel_state = vector {\n");
    model_label_print (file, mo, "\t", ",", ";");
  }
  fprintf (file, "\t};\n};\n\n");
}                               /* model_print */

/*============================================================================*/

void model_direct_print (FILE * file, model_direct * mo_d, int multip)
{
  int i, j;
  for (i = 0; i < multip; i++) {
    fprintf (file, "HMM = {\n\tM = %d;\n\tN = %d;\n", mo_d->M, mo_d->N);
    fprintf (file, "\tprior = %.3f;\n", mo_d->prior);
    fprintf (file, "\tA = matrix {\n");
    matrix_d_print (file, mo_d->A, mo_d->N, mo_d->N, "\t", ",", ";");
    fprintf (file, "\t};\n\tB = matrix {\n");
    matrix_d_print (file, mo_d->B, mo_d->N, mo_d->M, "\t", ",", ";");
    fprintf (file, "\t};\n\tPi = vector {\n");
    fprintf (file, "\t%.4f", mo_d->Pi[0]);
    for (j = 1; j < mo_d->N; j++)
      fprintf (file, ", %.4f", mo_d->Pi[j]);
    fprintf (file, ";\n\t};\n");
    fprintf (file, "\tfix_state = vector {\n");
    fprintf (file, "\t%d", mo_d->fix_state[0]);
    for (j = 1; j < mo_d->N; j++)
      fprintf (file, ", %d", mo_d->fix_state[j]);
    fprintf (file, ";\n\t};\n");
    fprintf (file, "};\n\n");
  }
}                               /* model_direct_print */

/*============================================================================*/

void model_direct_clean (model_direct * mo_d, hmm_check_t * check)
{
#define CUR_PROC "model_direct_clean"
  int i;
  if (!mo_d)
    return;
  mo_d->M = mo_d->N = 0;
  mo_d->prior = -1;
  if (mo_d->A) {
    for (i = 0; i < check->r_a; i++)
      m_free (mo_d->A[i]);
    m_free (mo_d->A);
  }
  if (mo_d->B) {
    for (i = 0; i < check->r_b; i++)
      m_free (mo_d->B[i]);
    m_free (mo_d->B);
  }
  if (mo_d->Pi){
    m_free (mo_d->Pi);
  }
  if (mo_d->fix_state){
    m_free (mo_d->fix_state);
  }
  
  mo_d->A = mo_d->B = NULL;
  mo_d->Pi = NULL;
  mo_d->fix_state = NULL;
#undef CUR_PROC
}                               /* model_direct_clean */

/*============================================================================*/

int model_direct_check_data (model_direct * mo_d, hmm_check_t * check)
{
#define CUR_PROC "model_direct_check_data"
  char *str;
  if (check->r_a != mo_d->N || check->c_a != mo_d->N) {
    str = mprintf (NULL, 0, "Incompatible dim. A (%d X %d) and N (%d)\n",
                   check->r_a, check->c_a, mo_d->N);
    mes_prot (str);
    m_free (str);
    return (-1);
  }
  if (check->r_b != mo_d->N || check->c_b != mo_d->M) {
    str =
      mprintf (NULL, 0, "Incompatible dim. B (%d X %d) and N X M (%d X %d)\n",
               check->r_b, check->c_b, mo_d->N, mo_d->M);
    mes_prot (str);
    m_free (str);
    return (-1);
  }
  if (check->len_pi != mo_d->N) {
    str = mprintf (NULL, 0, "Incompatible dim. Pi (%d) and N (%d)\n",
                   check->len_pi, mo_d->N);
    mes_prot (str);
    m_free (str);
    return (-1);
  }
  if (check->len_fix != mo_d->N) {
    str = mprintf (NULL, 0, "Incompatible dim. fix_state (%d) and N (%d)\n",
                   check->len_fix, mo_d->N);
    mes_prot (str);
    m_free (str);
    return (-1);
  }

  return 0;
#undef CUR_PROC
}                               /* model_direct_check_data */



/*============================================================================*/
/* XXX symmetric not implemented yet */
double model_prob_distance (model * m0, model * m, int maxT, int symmetric,
                            int verbose)
{
#define CUR_PROC "model_prob_distance"

#define STEPS 40

  double p0, p;
  double d = 0.0;
  double *d1;
  sequence_t *seq0 = NULL;
  sequence_t *tmp = NULL;
  model *mo1, *mo2;
  int i, t, a, k;
  int true_len;
  int true_number;
  int left_to_right = 0;
  long total, index;
  int step_width = 0;
  int steps = 1;

  /* printf("***  model_prob_distance:\n"); */

  if (verbose) {                /* If we are doing it verbosely we want to have 40 steps */
    step_width = maxT / 40;
    steps = STEPS;
  }
  else                          /* else just one */
    step_width = maxT;

  ARRAY_CALLOC (d1, steps);

  mo1 = m0;
  mo2 = m;

  for (k = 0; k < 2; k++) {     /* Two passes for the symmetric case */

    /* seed = 0 -> no reseeding. Call  ghmm_rng_timeseed(RNG) externally */
    seq0 = model_generate_sequences (mo1, 0, maxT + 1, 1, maxT + 1);



    if (seq0 == NULL) {
      mes_prot (" generate_sequences failed !");
      goto STOP;
    }

    if (seq0->seq_len[0] < maxT) {      /* There is an absorbing state */

      /* NOTA BENE: Assumpting the model delivers an explicit end state, 
         the condition of a fix initial state is removed. */

      /* For now check that Pi puts all weight on state */
      /*
         t = 0;
         for (i = 0; i < mo1->N; i++) {
         if (mo1->s[i].pi > 0.001)
         t++;
         }    
         if (t > 1) {
         mes_prot("ERROR: No proper left-to-right model. Multiple start states");
         goto STOP;
         } */

      left_to_right = 1;
      total = seq0->seq_len[0];

      while (total <= maxT) {

        /* create a additional sequences at once */
        a = (maxT - total) / (total / seq0->seq_number) + 1;
        /* printf("total=%d generating %d", total, a); */
        tmp = model_generate_sequences (mo1, 0, 0, a, a);
        if (tmp == NULL) {
          mes_prot (" generate_sequences failed !");
          goto STOP;
        }
        sequence_free (&tmp);
        sequence_add (seq0, tmp);

        total = 0;
        for (i = 0; i < seq0->seq_number; i++)
          total += seq0->seq_len[i];
      }
    }

    if (left_to_right) {

      for (t = step_width, i = 0; t <= maxT; t += step_width, i++) {

        index = 0;
        total = seq0->seq_len[0];

        /* Determine how many sequences we need to get a total of t
           and adjust length of last sequence to obtain total of 
           exactly t */

        while (total < t) {
          index++;
          total += seq0->seq_len[index];
        }

        true_len = seq0->seq_len[index];
        true_number = seq0->seq_number;

        if ((total - t) > 0)
          seq0->seq_len[index] = total - t;
        seq0->seq_number = index;

        p0 = model_likelihood (mo1, seq0);
        if (p0 == +1 || p0 == -1) {     /* error! */
          mes_prot ("problem: model_likelihood failed !");
          goto STOP;
        }
        p = model_likelihood (mo2, seq0);
        if (p == +1 || p == -1) {       /* what shall we do now? */
          mes_prot ("problem: model_likelihood failed !");
          goto STOP;
        }

        d = 1.0 / t * (p0 - p);

        if (symmetric) {
          if (k == 0)
            /* save d */
            d1[i] = d;
          else {
            /* calculate d */
            d = 0.5 * (d1[i] + d);
          }
        }

        if (verbose && (!symmetric || k == 1))
          printf ("%d\t%f\t%f\t%f\n", t, p0, p, d);

        seq0->seq_len[index] = true_len;
        seq0->seq_number = true_number;
      }
    }

    else {

      true_len = seq0->seq_len[0];

      for (t = step_width, i = 0; t <= maxT; t += step_width, i++) {
        seq0->seq_len[0] = t;

        p0 = model_likelihood (mo1, seq0);
        /* printf("   P(O|m1) = %f\n",p0); */
        if (p0 == +1) {         /* error! */
          mes_prot ("seq0 can't be build from mo1!");
          goto STOP;
        }
        p = model_likelihood (mo2, seq0);
        /* printf("   P(O|m2) = %f\n",p); */
        if (p == +1) {          /* what shall we do now? */
          mes_prot ("problem: seq0 can't be build from mo2!");
          goto STOP;
        }

        d = (1.0 / t) * (p0 - p);

        if (symmetric) {
          if (k == 0)
            /* save d */
            d1[i] = d;
          else {
            /* calculate d */
            d = 0.5 * (d1[i] + d);
          }
        }

        if (verbose && (!symmetric || k == 1))
          printf ("%d\t%f\t%f\t%f\n", t, p0, p, d);

      }
      seq0->seq_len[0] = true_len;
    }

    if (symmetric) {
      sequence_free (&seq0);
      mo1 = m;
      mo2 = m0;
    }
    else
      break;

  }                             /* k = 1,2 */

  sequence_free (&seq0);
  free (d1);
  return d;

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_free (&seq0);
  free (d1);
  return (0.0);
#undef CUR_PROC
}


/*============================================================================*/

void state_clean (state * my_state)
{
#define CUR_PROC "state_clean"

  if (my_state->b){
    m_free (my_state->b);
  }
  if (my_state->out_id){
    m_free (my_state->out_id);
  }
  if (my_state->in_id){
    m_free (my_state->in_id);
  }
  if (my_state->out_a){
    m_free (my_state->out_a);
  }
  if (my_state->in_a){
    m_free (my_state->in_a);
  }
  my_state->pi = 0;
  my_state->b = NULL;
  my_state->out_id = NULL;
  my_state->in_id = NULL;
  my_state->out_a = NULL;
  my_state->in_a = NULL;
  my_state->out_states = 0;
  my_state->in_states = 0;
  my_state->fix = 0;

#undef CUR_PROC
}                               /* state_clean */

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


 /*==========================Labeled HMMS ================================*/

sequence_t *model_label_generate_sequences (model * mo, int seed,
                                            int global_len, long seq_number,
                                            int Tmax)
{
# define CUR_PROC "model_label_generate_sequences"

  /* An end state is characterized by not having an output probabiliy. */

  sequence_t *sq = NULL;
  int i, j, m, transition_impossible, j_id;
  double p, sum;
  int len = global_len;
  /*int silent_len = 0; */
  int n = 0;
  int pos = 0;
  int label_index = 0;

  sq = sequence_calloc (seq_number);

  if (!sq) {
    mes_proc ();
    goto STOP;
  }

  /* allocating additional fields for the labels in the sequence_t struct */
  ARRAY_CALLOC (sq->state_labels, seq_number);
  ARRAY_CALLOC (sq->state_labels_len, seq_number);

  if (len <= 0)
    /* A specific length of the sequences isn't given. As a model should have
       an end state, the konstant MAX_SEQ_LEN is used. */
    len = (int) MAX_SEQ_LEN;

  if (seed > 0) {
    GHMM_RNG_SET (RNG, seed);
  }

  /* initialize the emission history */
  mo->emission_history = 0;

  while (n < seq_number) {
    /* printf("sequenz n = %d\n",n); */

    ARRAY_CALLOC (sq->seq[n], len);
    
    if (mo->model_type & kSilentStates) {
      
      /* printf("Model has silent states.\n"); */
      /* for silent models we have to allocate for the maximal possible number of lables */
      ARRAY_CALLOC (sq->state_labels[n], len * mo->N);
    }
    else {
      /* printf("Model has no silent states.\n"); */
      ARRAY_CALLOC (sq->state_labels[n], len);
    }
    label_index = 0;
    pos = 0;
    
    /* Get a random initial state i */
    p = GHMM_RNG_UNIFORM (RNG);
    sum = 0.0;
    for (i = 0; i < mo->N; i++) {
      sum += mo->s[i].pi;
      if (sum >= p)
        break;
    }

    if (!(mo->model_type & kHigherOrderEmissions) && 0 < mo->s[i].order) {
      fprintf (stderr,
               "ERROR: State %d has emission order %d, but it's initial probability is not 0.\n",
               i, mo->s[i].order);
      exit (1);
    }

    /* add label of fist state to the label list */
    sq->state_labels[n][label_index++] = mo->s[i].label;

    if (mo->model_type & kSilentStates && mo->silent[i]) {
      /* silent state: we do nothing, no output */
      /* printf("first state %d silent\n",i);
         pos = 0; 
         silent_len = silent_len + 1; */
    }
    else {
      /* first state emits */
      /* printf("first state %d not silent\n",i); */

      /* Get a random initial output m */
      m = get_random_output (mo, i, pos);
      update_emission_history (mo, m);
      sq->seq[n][pos++] = m;
    }

    /* check whether sequence was completed by inital state */
    if (pos >= len) {

      /* printf("assinging length %d to sequence %d\n", pos, n);
         printf("sequence complete...\n"); */
      sq->seq_len[n] = pos;
      sq->state_labels_len[n] = label_index;

      /* printf("1: seq %d -> %d labels\n",n,sq->state_labels_len[n]); */

      if (mo->model_type & kSilentStates) {
        printf ("reallocating\n");
        ARRAY_REALLOC (sq->state_labels[n], sq->state_labels_len[n]);
      }
      n++;
      continue;
    }
    while (pos < len) {
      if (pos < mo->maxorder) {
        /* if maxorder is not yet reached, we should only go to states with order < pos */
        transition_impossible = 1;
        for (j = 0; j < mo->s[i].out_states; j++) {
          j_id = mo->s[i].out_id[j];
          if ((mo->s[j_id].order) < (pos)) {
            transition_impossible = 0;
            break;
          }
        }
        if (1 == transition_impossible) {
          fprintf (stderr,
                   "No possible transition from state %d, due to too high order of all successor states. Position: %d",
                   i, pos);
          exit (1);
        }
        do {
          /* Get a random state i */
          p = GHMM_RNG_UNIFORM (RNG);
          sum = 0.0;
          for (j = 0; j < mo->s[i].out_states; j++) {
            j_id = mo->s[i].out_id[j];
            sum += mo->s[i].out_a[j];
            if (sum >= p)
              break;
          }
        } while (mo->s[j_id].order >= pos);
      }
      else {
        /* Get a random state i */
        p = GHMM_RNG_UNIFORM (RNG);
        sum = 0.0;
        for (j = 0; j < mo->s[i].out_states; j++) {
          sum += mo->s[i].out_a[j];
          if (sum >= p)
            break;
        }
      }
      i = mo->s[i].out_id[j];

      /* add label of state to the label list */
      sq->state_labels[n][label_index] = mo->s[i].label;
      label_index++;

      /* printf("state %d selected (i: %d, j: %d) at position %d\n",mo->s[i].out_id[j],i,j,pos); */

      if (sum == 0.0) {
	/* printf("final state reached - aborting\n"); */
	sq->seq_len[n++] = pos;
	break;
      }

      if (mo->model_type & kSilentStates && mo->silent[i]) {    /* Got a silent state i */
        /* printf("silent state \n");
           silent_len += 1;
           if (silent_len >= Tmax) {
           printf("%d silent states reached -> silent circle - aborting...\n",silent_len);
           sq->seq_len[n++] = pos; 
           break;
           } */
      }
      else {
        /* Get a random output m from state i */
        m = get_random_output (mo, i, pos);
        update_emission_history (mo, m);
        sq->seq[n][pos++] = m;
      }
      if (pos == len) {
        sq->state_labels_len[n] = label_index;

        /* printf("2: seq %d -> %d labels\n",n,sq->state_labels_len[n]); */

        if (mo->model_type & kSilentStates) {
          printf ("reallocating\n");
          ARRAY_REALLOC (sq->state_labels[n], sq->state_labels_len[n]);
        }

        sq->seq_len[n++] = pos;
      }
    }                           /* while (pos < len) */
  }                             /* while( n < seq_number ) */

  return (sq);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_free (&sq);
  return (NULL);
# undef CUR_PROC
}                               /* data */


/*----------------------------------------------------------------------------*/
/** gets correct index for emission array b of state j */
 int get_emission_index (model * mo, int j, int obs, int t)
{
  int retval;
  if (!(mo->model_type & kHigherOrderEmissions))
    return (obs);
  if (mo->s[j].order > t)
    retval = -1;
  else
    retval = (mo->emission_history * mo->M) % mo->pow_lookup[mo->s[j].order + 1] + obs;

  return retval;
}

/** updates emission history */
 void update_emission_history (model * mo, int obs)
{
  /* left-shift the history, truncate to history length and
     add the last observation */
  if (mo->model_type & kHigherOrderEmissions)
    mo->emission_history =
      (mo->emission_history * mo->M) % mo->pow_lookup[mo->maxorder] + obs;
}

/** updates emission history for backward algorithm*/
 void update_emission_history_front (model * mo, int obs)
{
  /*removes the most significant position (right-shift) and add the last seen
     observation (left-shifted with the length of history) */
  if (mo->model_type & kHigherOrderEmissions)
    mo->emission_history =
      (mo->pow_lookup[mo->maxorder - 1] * obs) + mo->emission_history / mo->M;
}


/*----------------------------------------------------------------------------*/
/* Scales the output and transitions probs of all states in a given model */
int model_normalize (model * mo)
{
#define CUR_PROC "model_normalize"

  int i, j, m, j_id, i_id=0, res=0;
  int size = 1;
  char *str;

  for (i = 0; i < mo->N; i++) {

    /* check model_type before using state order */
    if (mo->model_type & kHigherOrderEmissions)
      size = model_ipow (mo, mo->M, mo->s[i].order);

    /* normalize transition probabilities */
    if (vector_normalize (mo->s[i].out_a, mo->s[i].out_states) == -1) {
      res = -1;
    }
    /* for every outgoing probability update the corrosponding incoming probability */
    for (j = 0; j < mo->s[i].out_states; j++) {
      j_id = mo->s[i].out_id[j];
      for (m = 0; m < mo->s[j_id].in_states; m++) {
        if (i == mo->s[j_id].in_id[m]) {
          i_id = m;
          break;
        }
      }
      if (i_id == mo->s[j_id].in_states) {
        str = mprintf (NULL, 0, "Outgoing transition from state %d to \
           state %d has no corresponding incoming transition.\n", i, j_id);
        mes_prot (str);
        return -1;
      }
      mo->s[j_id].in_a[i_id] = mo->s[i].out_a[j];
    }
    /* normalize emission probabilities */
    for (m = 0; m < size; m++) {
      if (vector_normalize (&(mo->s[i].b[m * mo->M]), mo->M) == -1) {
        res = -1;
      }
    }
  }

  return res;
#undef CUR_PROC
}                               /* model_normalize */


/*----------------------------------------------------------------------------*/
int model_add_noise (model * mo, double level, int seed)
{
#define CUR_PROC "model_add_noise_A"

  int h, i, j, hist;
  int size = 1;

  if (level > 1.0)
    level = 1.0;

  for (i = 0; i < mo->N; i++) {
    for (j = 0; j < mo->s[i].out_states; j++)
      /* add noise only to out_a, in_a is updated on normalisation */
      mo->s[i].out_a[j] *= (1 - level) + (GHMM_RNG_UNIFORM (RNG) * 2 * level);

    if (mo->model_type & kHigherOrderEmissions)
      size = model_ipow (mo, mo->M, mo->s[i].order);
    for (hist = 0; hist < size; hist++)
      for (h = hist * mo->M; h < hist * mo->M + mo->M; h++)
        mo->s[i].b[h] *= (1 - level) + (GHMM_RNG_UNIFORM (RNG) * 2 * level);
  }

  return model_normalize (mo);

#undef CUR_PROC
}

/*----------------------------------------------------------------------------*/
int model_add_transition (state * s, int start, int dest, double prob)
{
#define CUR_PROC "model_add_transition"

  int i;

  /* resize the arrays */
  ARRAY_REALLOC (s[dest].in_id, s[dest].in_states + 1);
  ARRAY_REALLOC (s[dest].in_a, s[dest].in_states + 1);
  ARRAY_REALLOC (s[start].out_id, s[start].out_states + 1);
  ARRAY_REALLOC (s[start].out_a, s[start].out_states + 1);

  s[dest].in_states += 1;
  s[start].out_states += 1;

  /* search the right place to insert while moving greater entrys one field back */
  for (i = s[start].out_states - 1; i >= 0; i--) {
    if (i == 0 || dest > s[start].out_id[i - 1]) {
      s[start].out_id[i] = dest;
      s[start].out_a[i] = prob;
      break;
    }
    else {
      s[start].out_id[i] = s[start].out_id[i - 1];
      s[start].out_a[i] = s[start].out_a[i - 1];
    }
  }

  /* search the right place to insert while moving greater entrys one field back */
  for (i = s[dest].in_states - 1; i >= 0; i--)
    if (i == 0 || start > s[dest].in_id[i - 1]) {
      s[dest].in_id[i] = start;
      s[dest].in_a[i] = prob;
      break;
    }
    else {
      s[dest].in_id[i] = s[dest].in_id[i - 1];
      s[dest].in_a[i] = s[dest].in_a[i - 1];
    }

  return 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return -1;
#undef CUR_PROC
}

/*----------------------------------------------------------------------------*/
int model_del_transition (state * s, int start, int dest)
{
#define CUR_PROC "model_del_transition"

  int i, j;

  /* search ... */
  for (j = 0; dest != s[start].out_id[j]; j++)
    if (j == s[start].out_states) {
      mes_prot ("No such transition");
      return -1;
    }
  /* ... and replace outgoing */
  for (i = j + 1; i < s[start].out_states; i++) {
    s[start].out_id[i - 1] = s[start].out_id[i];
    s[start].out_a[i - 1] = s[start].out_a[i];
  }

  /* search ... */
  for (j = 0; start != s[dest].in_id[j]; j++)
    if (j == s[dest].in_states) {
      mes_prot ("No such transition");
      return -1;
    }
  /* ... and replace incoming */
  for (i = j + 1; i < s[dest].in_states; i++) {
    s[dest].in_id[i - 1] = s[dest].in_id[i];
    s[dest].in_a[i - 1] = s[dest].in_a[i];
  }

  /* reset number */
  s[dest].in_states -= 1;
  s[start].out_states -= 1;

  /* free memory */
  ARRAY_REALLOC (s[dest].in_id, s[dest].in_states);
  ARRAY_REALLOC (s[dest].in_a, s[dest].in_states);
  ARRAY_REALLOC (s[start].out_id, s[start].out_states);
  ARRAY_REALLOC (s[start].out_a, s[start].out_states);

  return 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return -1;
#undef CUR_PROC
}

/*----------------------------------------------------------------------------*/
/** 
   Allocates a new background_distributions struct and assigs the arguments to
   the respective fields. Note: The arguments need allocation outside of this
   function.
   
   @return     :               0 on success, -1 on error
   @param mo   :               one model
   @param cur  :               a id of a state
   @param times:               number of times the state cur is at least evaluated
*/
int model_apply_duration (model * mo, int cur, int times)
{
#define CUR_PROC "model_apply_duration"

  int i, j, last, size, failed=0;

  if (mo->model_type & kSilentStates) {
    mes_prot ("Sorry, apply_duration doesn't support silent states yet\n");
    return -1;
  }

  last = mo->N;
  mo->N += times - 1;

  ARRAY_REALLOC (mo->s, mo->N);
  if (mo->model_type & kSilentStates) {
    ARRAY_REALLOC (mo->silent, mo->N);
    ARRAY_REALLOC (mo->topo_order, mo->N);
  }
  if (mo->model_type & kTiedEmissions)
    ARRAY_REALLOC (mo->tied_to, mo->N);
  if (mo->model_type & kHasBackgroundDistributions)
    ARRAY_REALLOC (mo->background_id, mo->N);

  size = model_ipow (mo, mo->M, mo->s[cur].order + 1);
  for (i = last; i < mo->N; i++) {
    /* set the new state */
    mo->s[i].pi = 0.0;
    mo->s[i].order = mo->s[cur].order;
    mo->s[i].fix = mo->s[cur].fix;
    mo->s[i].label = mo->s[cur].label;
    mo->s[i].in_a = NULL;
    mo->s[i].in_id = NULL;
    mo->s[i].in_states = 0;
    mo->s[i].out_a = NULL;
    mo->s[i].out_id = NULL;
    mo->s[i].out_states = 0;

    ARRAY_MALLOC (mo->s[i].b, size);
    for (j = 0; j < size; j++)
      mo->s[i].b[j] = mo->s[cur].b[j];

    if (mo->model_type & kSilentStates) {
      mo->silent[i] = mo->silent[cur];
      /* XXX what to do with topo_order
         mo->topo_order[i] = ????????????; */
    }
    if (mo->model_type & kTiedEmissions)
      /* XXX is there a clean solution for tied states?
         what if the current state is a tie group leader?
         the last added state should probably become
         the new tie group leader */
      mo->tied_to[i] = kUntied;
    if (mo->model_type & kHasBackgroundDistributions)
      mo->background_id[i] = mo->background_id[cur];
  }

  /* move the outgoing transitions to the last state */
  while (mo->s[cur].out_states > 0) {
    if (mo->s[cur].out_id[0] == cur) {
      model_add_transition (mo->s, mo->N - 1, mo->N - 1, mo->s[cur].out_a[0]);
      model_del_transition (mo->s, cur, mo->s[cur].out_id[0]);
    }
    else {
      model_add_transition (mo->s, mo->N - 1, mo->s[cur].out_id[0],
                            mo->s[cur].out_a[0]);
      model_del_transition (mo->s, cur, mo->s[cur].out_id[0]);
    }
  }

  /* set the linear transitions through all added states */
  model_add_transition (mo->s, cur, last, 1.0);
  for (i = last + 1; i < mo->N; i++) {
    model_add_transition (mo->s, i - 1, i, 1.0);
  }

  if (model_normalize (mo))
    goto STOP;

  return 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  /* Fail hard if these realloc fail. They shouldn't because we have the memory
     and try only to clean up! */
  if (failed++)
    exit (1);
  
  ARRAY_REALLOC (mo->s, last);
  ARRAY_REALLOC (mo->tied_to, last);
  ARRAY_REALLOC (mo->background_id, last);

  mo->N = last;
  return -1;
#undef CUR_PROC
}

/*----------------------------------------------------------------------------*/
/** 
   Allocates a new background_distributions struct and assigs the arguments to
   the respective fields. Note: The arguments need allocation outside of this
   function.
   
   @return    :               new pointer to a background_distributions struct
   @param n   :               number of distributions
   @param order:              orders of the distribtions
   @param B:                  matrix of distribution parameters
*/
background_distributions *model_alloc_background_distributions (int n, int m,
                                                                int *orders,
                                                                double **B)
{
#define CUR_PROC "model_alloc_background_distributions"
  background_distributions *ptbackground;

  ARRAY_CALLOC (ptbackground, 1);

  ptbackground->n = n;
  ptbackground->m = m;
  if (orders != NULL && B != NULL) {
    ptbackground->order = orders;
    ptbackground->b = B;
  }
  else {
    m_free (ptbackground);
    return NULL;
  }
  return ptbackground;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return NULL;
#undef CUR_PROC
}

/*----------------------------------------------------------------------------*/

background_distributions
  *model_copy_background_distributions (background_distributions * bg)
{
#define CUR_PROC "model_copy_background_distributions"
  int i, j, b_i_len;
  int *new_order;
  double **new_b;

  ARRAY_MALLOC (new_order, bg->n);
  ARRAY_CALLOC (new_b, bg->n);

  for (i = 0; i < bg->n; i++) {
    new_order[i] = bg->order[i];
    b_i_len = pow (bg->m, bg->order[i] + 1);
    ARRAY_CALLOC (new_b[i], b_i_len);
    for (j = 0; j < b_i_len; j++) {
      new_b[i][j] = bg->b[i][j];
    }
  }

  return model_alloc_background_distributions (bg->n, bg->m, new_order,
                                               new_b);

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */

  return NULL;

#undef CUR_PROC
}


/*----------------------------------------------------------------------------*/
int model_apply_background (model * mo, double *background_weight)
{
# define CUR_PROC "model_apply_background"

  int i, j, size;

  if (!mo->model_type && kHasBackgroundDistributions) {
    mes_prot ("Error: No background distributions");
    return -1;
  }

  for (i = 0; i < mo->N; i++) {
    if (mo->background_id[i] != kNoBackgroundDistribution) {
      if (mo->s[i].order != mo->bp->order[mo->background_id[i]]) {
        mes_prot ("Error: State and background order do not match\n");
        return -1;
      }

      /* XXX Cache in background_distributions */
      size = model_ipow (mo, mo->M, mo->s[i].order + 1);
      for (j = 0; j < size; j++)
        mo->s[i].b[j] = (1.0 - background_weight[i]) * mo->s[i].b[j]
          + background_weight[i] * mo->bp->b[mo->background_id[i]][j];
    }
  }

  return 0;
#undef CUR_PROC
}                               /* model_apply_background */


/*----------------------------------------------------------------------------*/
int model_get_uniform_background (model * mo, sequence_t * sq)
{
# define CUR_PROC "get_background"

  int h, i, j, m, t, n=0;
  int e_index, size;
  double sum=0.0;

  if (!(mo->model_type & kHasBackgroundDistributions)) {
    mes_prot ("Error: Model has no background distribution");
    return -1;
  }

  mo->bp = NULL;
  ARRAY_MALLOC (mo->background_id, mo->N);

  /* create a background distribution for each state */
  for (i = 0; i < mo->N; i++) {
    mo->background_id[i] = mo->s[i].order;
  }

  /* allocate */
  ARRAY_CALLOC (mo->bp, 1);
  ARRAY_CALLOC (mo->bp->order, mo->maxorder);

  /* set number of distributions */
  mo->bp->n = mo->maxorder;

  /* set br->order */
  for (i = 0; i < mo->N; i++)
    if (mo->background_id[i] != kNoBackgroundDistribution)
      mo->bp->order[mo->background_id[i]] = mo->s[i].order;

  /* allocate and initialize br->b with zeros */
  ARRAY_CALLOC (mo->bp->b, mo->bp->n);

  for (i = 0; i < mo->bp->n; i++)
    ARRAY_MALLOC (mo->bp->b[i], model_ipow (mo, mo->M, mo->bp->order[i] + 1));

  for (i = 0; i < mo->bp->n; i++) {

    /* find a state with the current order */
    for (j = 0; j < mo->N; j++)
      if (mo->bp->order[i] == mo->s[j].order)
        break;

    /* initialize with ones as psoudocounts */
    size = model_ipow (mo, mo->M, mo->bp->order[n] + 1);
    for (m = 0; m < size; m++)
      mo->bp->b[i][m] = 1.0;

    for (n = 0; n < sq->seq_number; n++) {

      for (t = 0; t < mo->bp->order[i]; t++)
        update_emission_history (mo, sq->seq[n][t]);

      for (t = mo->bp->order[i]; t < sq->seq_len[n]; t++) {

        e_index = get_emission_index (mo, j, sq->seq[n][t], t);
        if (-1 != e_index)
          mo->bp->b[i][e_index]++;
      }
    }

    /* normalise */
    size = model_ipow (mo, mo->M, mo->bp->order[n]);
    for (h = 0; h < size; h += mo->M) {
      for (m = h; m < h + mo->M; m++)
        sum += mo->bp->b[i][m];
      for (m = h; m < h + mo->M; m++)
        mo->bp->b[i][m] /= sum;
    }

  }

  return 0;

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */


  return -1;
# undef CUR_PROC
}                               /* end get_background */


double model_distance(const model * mo, const model * m2) {
#define CUR_PROC "model_distances"

  int i, j, number=0;
  double tmp, distance=0.0;

/*   if (!model_check_compatibility(mo, m2)) */
/*     exit(1); */
/*   if (!model_check(mo)) */
/*     exit(1); */
/*   if (!model_check(m2)) */
/*     exit(1); */


  /* PI */
  for (i=0; i<mo->N; ++i) {
    tmp = mo->s[i].pi - m2->s[i].pi;
    distance += tmp*tmp;
    ++number;
  }
  for (i=0; i<mo->N; ++i) {
    /* A */
    for (j=0; j<mo->s[i].out_states; ++j) {
      tmp = mo->s[i].out_a[j] - m2->s[i].out_a[j];
      distance += tmp*tmp;
      ++number;
    }
    /* B */
    for (j=0; j<model_ipow(mo, mo->M, mo->s[i].order+1); ++j) {
      tmp = mo->s[i].b[j] - m2->s[i].b[j];
      distance += tmp*tmp;
      ++number;
    } 
  }

  return (distance/number);
#undef CUR_PROC
}

/*============================================================================*/

/*===================== E n d   o f  f i l e  "model.c"       ===============*/
