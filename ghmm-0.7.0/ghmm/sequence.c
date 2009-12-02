/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/sequence.c
*       Authors:  Bernd Wichern, Andrea Weisse, Utz J. Pape, Benjamin Georgi
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

#include "mprintf.h"
#include "mes.h"
#include "sequence.h"
#include <math.h>
#include <float.h>
#include "mprintf.h"
#include "mes.h"
#include "matrix.h"
#include "vector.h"
#include "const.h"
#include "model.h"
#include "foba.h"
#include "sfoba.h"
#include "vector.h"
#include "rng.h"
#include "string.h"
#include <ghmm/internal.h>

/*============================================================================*/
sequence_t **sequence_read (const char *filename, int *sq_number)
{
#define CUR_PROC "sequence_read"
  int i;
  sequence_t **sequence = NULL;
  scanner_t *s = NULL;
  *sq_number = 0;
  s = scanner_alloc (filename);
  if (!s) {
    mes_proc ();
    goto STOP;
  }

  while (!s->err && !s->eof && s->c - '}') {
    scanner_get_name (s);
    scanner_consume (s, '=');
    if (s->err)
      goto STOP;
    /* sequence file */
    if (!strcmp (s->id, "SEQ")) {
      (*sq_number)++;
      /* more mem */
      ARRAY_REALLOC (sequence, *sq_number);
      sequence[*sq_number - 1] = sequence_read_alloc (s);
      if (!sequence[*sq_number - 1]) {
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
  }

  scanner_free (&s);
  return sequence;


STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  scanner_free (&s);
  for (i = 0; i < *sq_number; i++)
    sequence_free (&(sequence[i]));
  m_free (sequence);
  *sq_number = 0;
  return NULL;
#undef CUR_PROC
}

/*============================================================================*/
sequence_t *sequence_read_alloc (scanner_t * s)
{
#define CUR_PROC "sequence_read_alloc"
  int symbols = 0, lexWord = 0;
  sequence_t *sq = NULL;
  int seq_len_lex = 0;
  ARRAY_CALLOC (sq, 1);
  scanner_consume (s, '{');
  if (s->err)
    goto STOP;
  while (!s->err && !s->eof && s->c - '}') {
    scanner_get_name (s);
    scanner_consume (s, '=');
    if (s->err)
      goto STOP;
    /* array of sequences to read */
    if (!strcmp (s->id, "O")) {
      scanner_consume (s, '{');
      if (s->err)
        goto STOP;
      sq->seq_number = 0;
      sq->total_w = 0.0;
      while (!s->eof && !s->err && s->c - '}') {
        /* another sequence --> realloc */
        ARRAY_REALLOC (sq->seq, sq->seq_number + 1);
        ARRAY_REALLOC (sq->seq_len, sq->seq_number + 1);
        ARRAY_REALLOC (sq->seq_label, sq->seq_number + 1);
        ARRAY_REALLOC (sq->seq_id, sq->seq_number + 1);
        ARRAY_REALLOC (sq->seq_w, sq->seq_number + 1);
        /* Label and ID */
        /* default */
        sq->seq_label[sq->seq_number] = -1;
        sq->seq_id[sq->seq_number] = -1.0;
        sq->seq_w[sq->seq_number] = 1;
        while (s->c == '<' || s->c == '(' || s->c == '|') {
          if (s->c == '<') {
            scanner_consume (s, '<');
            if (s->err)
              goto STOP;
            sq->seq_label[sq->seq_number] = scanner_get_int (s);
            if (s->err)
              goto STOP;
            scanner_consume (s, '>');
            if (s->err)
              goto STOP;
          }
          if (s->c == '(') {
            scanner_consume (s, '(');
            if (s->err)
              goto STOP;
            sq->seq_id[sq->seq_number] = scanner_get_edouble (s);
            if (s->err)
              goto STOP;
            scanner_consume (s, ')');
            if (s->err)
              goto STOP;
          }
          if (s->c == '|') {
            scanner_consume (s, '|');
            if (s->err)
              goto STOP;
            sq->seq_w[sq->seq_number] = (double) scanner_get_int (s);
            if (sq->seq_w[sq->seq_number] <= 0) {
              scanner_error (s, "sequence weight not positiv\n");
              goto STOP;
            }
            if (s->err)
              goto STOP;
            scanner_consume (s, '|');
            if (s->err)
              goto STOP;
          }
        }

        sq->seq[sq->seq_number] =
          scanner_get_int_array (s, sq->seq_len + sq->seq_number);
        if (sq->seq_len[sq->seq_number] > MAX_SEQ_LEN) {
          scanner_error (s, "sequence too long");
          goto STOP;
        }
        scanner_consume (s, ';');
        if (s->err)
          goto STOP;
        sq->total_w += sq->seq_w[sq->seq_number];
        sq->seq_number++;
      }                         /* while( !s->eof...) */
      if ((sq->seq_number == 0) || (sq->seq_number > MAX_SEQ_NUMBER)) {
        char *str = mprintf (NULL, 0,
                             "Number of sequences %ld exceeds possible range",
                             sq->seq_number);
        mes_prot (str);
        m_free (str);
        goto STOP;
      }
      scanner_consume (s, '}');
      if (s->err)
        goto STOP;
    }
    /* all possible seqs., sorted lexicographical */
    else if (!strcmp (s->id, "L")) {
      lexWord = 1;
      scanner_consume (s, '{');
      if (s->err)
        goto STOP;
      while (!s->err && !s->eof && s->c - '}') {
        scanner_get_name (s);
        scanner_consume (s, '=');
        if (s->err)
          goto STOP;
        if (!strcmp (s->id, "seq_len")) {
          seq_len_lex = scanner_get_int (s);
          if (s->err)
            goto STOP;
          if (seq_len_lex <= 0) {
            mes_prot ("Value for sequence length not allowed");
            goto STOP;
          }
        }
        else if (!strcmp (s->id, "symb")) {
          if (symbols < 0) {
            mes_prot ("Value for number of symbols not allowed");
            goto STOP;
          }
          symbols = scanner_get_int (s);
          if (s->err)
            goto STOP;
        }
        else {
          scanner_error (s, "unknown identifier");
          goto STOP;
        }
        scanner_consume (s, ';');
        if (s->err)
          goto STOP;
      }
      scanner_consume (s, '}');
      if ((seq_len_lex <= 0) || (symbols < 0)) {
        mes_prot
          ("Values for seq. length or number of symbols not spezified");
        goto STOP;
      }
      sq = sequence_lexWords (seq_len_lex, symbols);
      if (!sq)
        goto STOP;
    }                           /*if (!strcmp(s->id, "L")) */
    else {
      scanner_error (s, "unknown identifier");
      goto STOP;
    }
    scanner_consume (s, ';');
    if (s->err)
      goto STOP;
  }                             /* while(!s->err && !s->eof && s->c - '}') */
  scanner_consume (s, '}');
  if (s->err)
    goto STOP;
  return (sq);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_free (&sq);
  return (NULL);
#undef CUR_PROC
}                               /* sequence_read_alloc */

/*============================================================================*/

sequence_d_t **sequence_d_read (const char *filename, int *sqd_number)
{
#define CUR_PROC "sequence_d_read"
  int i;
  scanner_t *s = NULL;
  sequence_d_t **sequence = NULL;
  *sqd_number = 0;
  s = scanner_alloc (filename);
  if (!s) {
    mes_proc ();
    goto STOP;
  }

  while (!s->err && !s->eof && s->c - '}') {
    scanner_get_name (s);
    scanner_consume (s, '=');
    if (s->err)
      goto STOP;
    /* sequence file */
    if (!strcmp (s->id, "SEQD")) {
      (*sqd_number)++;
      /* more mem */
      ARRAY_REALLOC (sequence, *sqd_number);
      sequence[*sqd_number - 1] = sequence_d_read_alloc (s);
      if (!sequence[*sqd_number - 1]) {
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
  }
  scanner_free (&s);

  return sequence;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  scanner_free (&s);
  for (i = 0; i < *sqd_number; i++)
    sequence_d_free (&(sequence[i]));
  m_free (sequence);
  *sqd_number = 0;
  
  
  return NULL;
#undef CUR_PROC
}                               /* sequence_d_read */

/*============================================================================*/
sequence_d_t *sequence_d_read_alloc (scanner_t * s)
{
#define CUR_PROC "sequence_d_read_alloc"
  sequence_d_t *sqd = NULL;
  ARRAY_CALLOC (sqd, 1);
  scanner_consume (s, '{');
  if (s->err)
    goto STOP;
  while (!s->err && !s->eof && s->c - '}') {
    scanner_get_name (s);
    scanner_consume (s, '=');
    if (s->err)
      goto STOP;
    /* array of sequences to read */
    if (!strcmp (s->id, "O")) {
      scanner_consume (s, '{');
      if (s->err)
        goto STOP;
      sqd->seq_number = 0;
      sqd->total_w = 0.0;
      while (!s->eof && !s->err && s->c - '}') {
        /* another sequence --> realloc */
        ARRAY_REALLOC (sqd->seq, sqd->seq_number + 1);
        ARRAY_REALLOC (sqd->seq_len, sqd->seq_number + 1);
        ARRAY_REALLOC (sqd->seq_label, sqd->seq_number + 1);
        ARRAY_REALLOC (sqd->seq_id, sqd->seq_number + 1);
        ARRAY_REALLOC (sqd->seq_w, sqd->seq_number + 1);
        /* Label and ID and weight */
        /* default */
        sqd->seq_label[sqd->seq_number] = -1;
        sqd->seq_id[sqd->seq_number] = -1.0;
        sqd->seq_w[sqd->seq_number] = 1;
        while (s->c == '<' || s->c == '(' || s->c == '|') {
          if (s->c == '<') {
            scanner_consume (s, '<');
            if (s->err)
              goto STOP;
            sqd->seq_label[sqd->seq_number] = scanner_get_int (s);
            if (s->err)
              goto STOP;
            scanner_consume (s, '>');
            if (s->err)
              goto STOP;
          }
          if (s->c == '(') {
            scanner_consume (s, '(');
            if (s->err)
              goto STOP;
            sqd->seq_id[sqd->seq_number] = scanner_get_edouble (s);
            if (s->err)
              goto STOP;
            scanner_consume (s, ')');
            if (s->err)
              goto STOP;
          }
          if (s->c == '|') {
            scanner_consume (s, '|');
            if (s->err)
              goto STOP;
            sqd->seq_w[sqd->seq_number] = (double) scanner_get_int (s);
            if (sqd->seq_w[sqd->seq_number] < 0) {
              scanner_error (s, "negativ sequence weight\n");
              goto STOP;
            }
            if (s->err)
              goto STOP;
            scanner_consume (s, '|');
            if (s->err)
              goto STOP;
          }
        }
        sqd->seq[sqd->seq_number] =
          scanner_get_double_earray (s, sqd->seq_len + sqd->seq_number);
        if (sqd->seq_len[sqd->seq_number] > MAX_SEQ_LEN) {
          scanner_error (s, "sequence too long");
          goto STOP;
        }
        scanner_consume (s, ';');
        if (s->err)
          goto STOP;
        sqd->total_w += sqd->seq_w[sqd->seq_number];
        sqd->seq_number++;
      }                         /* while( !s->eof...) */
      if ((sqd->seq_number == 0) || (sqd->seq_number > MAX_SEQ_NUMBER)) {
        char *str = mprintf (NULL, 0,
                             "Number of sequences %ld exceeds possible range",
                             sqd->seq_number);
        mes_prot (str);
        m_free (str);
        goto STOP;
      }
      scanner_consume (s, '}');
      if (s->err)
        goto STOP;
    }
    else {
      scanner_error (s, "unknown identifier");
      goto STOP;
    }
    scanner_consume (s, ';');
    if (s->err)
      goto STOP;
  }                             /* while(!s->err && !s->eof && s->c - '}') */
  scanner_consume (s, '}');
  if (s->err)
    goto STOP;
  return (sqd);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_d_free (&sqd);
  return (NULL);
#undef CUR_PROC
}                               /* sequence_d_read_alloc */

/*============================================================================*/

/* Truncate Sequences in a given sqd_field; useful for Testing;
   returns truncated sqd_field; 
   trunc_ratio 0: no truncation
   trunc_ratio 1: truncation (mean truncation faktor = 0.5)
   trunc_ratio -1: 100 % truncation
*/

sequence_d_t **sequence_d_truncate (sequence_d_t ** sqd_in, int sqd_fields,
                                    double trunc_ratio, int seed)
{
#define CUR_PROC "sequence_d_truncate"
  sequence_d_t **sq;
  int i, j, trunc_len;
  /* Hack, use -1 for complete truncation */
  if ((0 > trunc_ratio || 1 < trunc_ratio) && trunc_ratio != -1) {
    mes_prot ("Error: trunc_ratio not valid\n");
    goto STOP;
  }
  ARRAY_CALLOC (sq, sqd_fields);

  ghmm_rng_init ();
  GHMM_RNG_SET (RNG, seed);

  for (i = 0; i < sqd_fields; i++) {
    sq[i] = sequence_d_calloc (sqd_in[i]->seq_number);
    sq[i]->total_w = sqd_in[i]->total_w;
    for (j = 0; j < sqd_in[i]->seq_number; j++) {
      ARRAY_CALLOC (sq[i]->seq[j], sqd_in[i]->seq_len[j]);
      /* length of truncated seq. */
      if (trunc_ratio == -1)
        trunc_len = 0;
      else
        trunc_len = (int) ceil ((sqd_in[i]->seq_len[j] *
                                 (1 - trunc_ratio * GHMM_RNG_UNIFORM (RNG))));
      sequence_d_copy (sq[i]->seq[j], sqd_in[i]->seq[j], trunc_len);
      ARRAY_REALLOC (sq[i]->seq[j], trunc_len);
      sq[i]->seq_len[j] = trunc_len;
      sq[i]->seq_label[j] = sqd_in[i]->seq_label[j];
      sq[i]->seq_id[j] = sqd_in[i]->seq_id[j];
      sq[i]->seq_w[j] = sqd_in[i]->seq_w[j];
    }
  }                             /* for all sqd_fields */

  return sq;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return NULL;
#undef CUR_PROC
}

/*============================================================================*/
sequence_d_t *sequence_d_calloc (long seq_number)
{
#define CUR_PROC "sequence_dcalloc"
  int i;
  /*printf("*** sequence_d_t *sequence_d_calloc, nr: %d\n",seq_number);*/
  sequence_d_t *sqd = NULL;
  if (seq_number > MAX_SEQ_NUMBER) {
    char *str = mprintf (NULL, 0,
                         "Number of sequences %ld exceeds possible range",
                         seq_number);
    mes_prot (str);
    m_free (str);
    goto STOP;
  }
  ARRAY_CALLOC (sqd, 1);
  ARRAY_CALLOC (sqd->seq, seq_number);
  ARRAY_CALLOC (sqd->seq_len, seq_number);
  ARRAY_CALLOC (sqd->seq_label, seq_number);
  ARRAY_CALLOC (sqd->seq_id, seq_number);
  ARRAY_CALLOC (sqd->seq_w, seq_number);
  sqd->seq_number = seq_number;

  sqd->total_w = 0.0;
  for (i = 0; i < seq_number; i++) {
    sqd->seq_label[i] = -1;
    sqd->seq_id[i] = -1.0;
    sqd->seq_w[i] = 1;
  }

  return sqd;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_d_free (&sqd);
  return NULL;
#undef CUR_PROC
}                               /* sequence_d_calloc */

/*============================================================================*/
sequence_t *sequence_calloc (long seq_number)
{
#define CUR_PROC "sequence_calloc"
  int i;
  sequence_t *sq = NULL;
  if (seq_number > MAX_SEQ_NUMBER) {
    char *str = mprintf (NULL, 0,
                         "Number of sequences %ld exceeds possible range",
                         seq_number);
    mes_prot (str);
    m_free (str);
    goto STOP;
  }
  ARRAY_CALLOC (sq, 1);
  ARRAY_CALLOC (sq->seq, seq_number);
  /*ARRAY_CALLOC (sq->states, seq_number);*/
  ARRAY_CALLOC (sq->seq_len, seq_number);
  ARRAY_CALLOC (sq->seq_label, seq_number);
  ARRAY_CALLOC (sq->seq_id, seq_number);
  ARRAY_CALLOC (sq->seq_w, seq_number);
  sq->seq_number = seq_number;
  for (i = 0; i < seq_number; i++) {
    sq->seq_label[i] = -1;
    sq->seq_id[i] = -1.0;
    sq->seq_w[i] = 1;
  }
  return sq;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_free (&sq);
  return NULL;
#undef CUR_PROC
}                               /* sequence_calloc */

/*============================================================================*/

sequence_d_t *sequence_d_get_singlesequence(sequence_d_t *sq, int index)
{
  sequence_d_t *res;
  res = sequence_d_calloc(1);
  
  res->seq[0] = sq->seq[index];
  res->seq_len[0] = sq->seq_len[index];
  res->seq_label[0] = sq->seq_label[index];
  res->seq_id[0] = sq->seq_id[index];
  res->seq_w[0] = sq->seq_w[index];
  res->total_w = res->seq_w[0];

  return res;
  
}

sequence_t *sequence_get_singlesequence(sequence_t *sq, int index)
{
#define CUR_PROC "sequence_get_singlesequence"
  sequence_t *res;
  res = sequence_calloc(1);
  if (!res) goto STOP;
  
  res->seq[0] = sq->seq[index];
  res->seq_len[0] = sq->seq_len[index];
  res->seq_label[0] = sq->seq_label[index];
  res->seq_id[0] = sq->seq_id[index];
  res->seq_w[0] = sq->seq_w[index];
  res->total_w = res->seq_w[0];

  if (sq->state_labels){
      ARRAY_CALLOC (res->state_labels, 1);
      ARRAY_CALLOC (res->state_labels_len, 1);
      res->state_labels[0] = sq->state_labels[index];
      res->state_labels_len[0] = sq->state_labels_len[index];
  }
  
  return res;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return NULL;
#undef CUR_PROC
}
/*XXX TEST: frees everything but the seq field */
int sequence_subseq_free (sequence_t ** sq)
{
# define CUR_PROC "sequence_subseq_free"
  /*int i,j;*/

  mes_check_ptr (sq, return (-1));
  if (!*sq)
    return (0);

  m_free ((*sq)->seq_len);
  m_free ((*sq)->seq_label);
  m_free ((*sq)->seq_id);
  m_free ((*sq)->seq_w);

  if ((*sq)->states) {
    matrix_i_free (&(*sq)->states, (*sq)->seq_number);
    /*m_free((*sq)->states); */
  }

  if ((*sq)->state_labels) {
    matrix_i_free (&(*sq)->state_labels, (*sq)->seq_number);
    m_free ((*sq)->state_labels_len);

    /*m_free((*sq)->states); */
  }

  m_free ((*sq)->seq);
  m_free (*sq);
  return 0;
# undef CUR_PROC
}                               /* sequence_subseq_free */


int sequence_d_subseq_free (sequence_d_t ** sqd)
{
# define CUR_PROC "sequence_d_subseq_free"
  mes_check_ptr (sqd, return (-1));


  /* sequence_d_print(stdout,*sqd,0);*/
  m_free ((*sqd)->seq);
  m_free ((*sqd)->seq_len);
  m_free ((*sqd)->seq_label);
  m_free ((*sqd)->seq_id);
  m_free ((*sqd)->seq_w);
  m_free (*sqd);
  return 0;
# undef CUR_PROC
}   /* sequence_d_subseq_free */



/*============================================================================*/
sequence_t *sequence_lexWords (int n, int M)
{
# define CUR_PROC "sequence_lexWords"

  sequence_t *sq = NULL;
  long seq_number, cnt = 0;
  int j = n - 1;
  int i;
  int *seq;
  if ((n < 0) || (M <= 0)) {
    mes_proc ();
    goto STOP;
  }
  seq_number = (long) pow ((double) M, (double) n);
  sq = sequence_calloc (seq_number);
  if (!sq) {
    mes_proc ();
    goto STOP;
  }
  for (i = 0; i < seq_number; i++) {
    ARRAY_CALLOC (sq->seq[i], n);
    sq->seq_len[i] = n;
    sq->seq_id[i] = i;
  }

  ARRAY_CALLOC (seq, n);
  while (!(j < 0)) {
    sequence_copy (sq->seq[cnt], seq, n);
    j = n - 1;
    while (seq[j] == M - 1) {
      seq[j] = 0;
      j--;
    }
    seq[j]++;
    cnt++;
  }
  m_free (seq);
  return sq;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_free (&sq);
  return NULL;
# undef CUR_PROC
}                               /* sequence_lewWords */

/*============================================================================*/
int sequence_max_symbol (sequence_t * sq)
{
  long i, j;
  int max_symb = -1;
  for (i = 0; i < sq->seq_number; i++)
    for (j = 0; j < sq->seq_len[i]; j++) {
      if (sq->seq[i][j] > max_symb)
        max_symb = sq->seq[i][j];
    }
  return max_symb;
}                               /* sequence_max_symbol */

/*============================================================================*/
void sequence_copy (int *target, int *source, int len)
{
  int i;
  for (i = 0; i < len; i++)
    target[i] = source[i];
}                               /* sequence_copy */


/*============================================================================*/
void sequence_d_copy (double *target, double *source, int len)
{
  int i;
  for (i = 0; i < len; i++)
    target[i] = source[i];
}                               /* sequence_copy */


/*============================================================================*/
int sequence_add (sequence_t * target, sequence_t * source)
{
#define CUR_PROC "sequence_add"

  int res = -1;
  int **old_seq = target->seq;
  /*int **old_seq_st    = target->states;*/
  int *old_seq_len = target->seq_len;
  long *old_seq_label = target->seq_label;
  double *old_seq_id = target->seq_id;
  double *old_seq_w = target->seq_w;
  long old_seq_number = target->seq_number;
  long i;

  target->seq_number = old_seq_number + source->seq_number;
  target->total_w += source->total_w;

  ARRAY_CALLOC (target->seq, target->seq_number);
  /*ARRAY_CALLOC (target->states, target->seq_number);*/
  ARRAY_CALLOC (target->seq_len, target->seq_number);
  ARRAY_CALLOC (target->seq_label, target->seq_number);
  ARRAY_CALLOC (target->seq_id, target->seq_number);
  ARRAY_CALLOC (target->seq_w, target->seq_number);

  for (i = 0; i < old_seq_number; i++) {
    target->seq[i] = old_seq[i];
    /*target->states[i] = old_seq_st[i];*/
    target->seq_len[i] = old_seq_len[i];
    target->seq_label[i] = old_seq_label[i];
    target->seq_id[i] = old_seq_id[i];
    target->seq_w[i] = old_seq_w[i];
  }

  for (i = 0; i < (target->seq_number - old_seq_number); i++) {
    ARRAY_CALLOC (target->seq[i + old_seq_number], source->seq_len[i]);

    sequence_copy (target->seq[i + old_seq_number], source->seq[i],
                   source->seq_len[i]);

    /*ARRAY_CALLOC (target->states[i+old_seq_number], source->seq_len[i]); */

    /* sequence_copy(target->states[i+old_seq_number], source->states[i], 
       source->seq_len[i]); */

    target->seq_len[i + old_seq_number] = source->seq_len[i];
    target->seq_label[i + old_seq_number] = source->seq_label[i];
    target->seq_id[i + old_seq_number] = source->seq_id[i];
    target->seq_w[i + old_seq_number] = source->seq_w[i];
  }


  m_free (old_seq);
  /*m_free(old_seq_st);*/
  m_free (old_seq_len);
  m_free (old_seq_label);
  m_free (old_seq_id);
  m_free (old_seq_w);
  res = 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return res;
#undef CUR_PROC
}


/*============================================================================*/
int sequence_d_add (sequence_d_t * target, sequence_d_t * source)
{
#define CUR_PROC "sequence_d_add"

  int res = -1;
  double **old_seq = target->seq;
  int *old_seq_len = target->seq_len;
  long *old_seq_label = target->seq_label;
  double *old_seq_id = target->seq_id;
  double *old_seq_w = target->seq_w;
  long old_seq_number = target->seq_number;
  long i;

  target->seq_number = old_seq_number + source->seq_number;
  target->total_w += source->total_w;

  ARRAY_CALLOC (target->seq, target->seq_number);
  ARRAY_CALLOC (target->seq_len, target->seq_number);
  ARRAY_CALLOC (target->seq_label, target->seq_number);
  ARRAY_CALLOC (target->seq_id, target->seq_number);
  ARRAY_CALLOC (target->seq_w, target->seq_number);

  for (i = 0; i < old_seq_number; i++) {
    target->seq[i] = old_seq[i];
    target->seq_len[i] = old_seq_len[i];
    target->seq_label[i] = old_seq_label[i];
    target->seq_id[i] = old_seq_id[i];
    target->seq_w[i] = old_seq_w[i];
  }

  for (i = 0; i < (target->seq_number - old_seq_number); i++) {
    ARRAY_CALLOC (target->seq[i + old_seq_number], source->seq_len[i]);

    sequence_d_copy (target->seq[i + old_seq_number], source->seq[i],
                     source->seq_len[i]);
    target->seq_len[i + old_seq_number] = source->seq_len[i];
    target->seq_label[i + old_seq_number] = source->seq_label[i];
    target->seq_id[i + old_seq_number] = source->seq_id[i];
    target->seq_w[i + old_seq_number] = source->seq_w[i];
  }

  m_free (old_seq);
  m_free (old_seq_len);
  m_free (old_seq_label);
  m_free (old_seq_id);
  m_free (old_seq_w);
  res = 0;

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return res;
#undef CUR_PROC
}

/*============================================================================*/
int sequence_check (sequence_t * sq, int max_symb)
{
#define CUR_PROC "sequence_check"
  int i, j;
  for (j = 0; j < sq->seq_number; j++) {
    for (i = 0; i < sq->seq_len[j]; i++) {
      if ((sq->seq[j][i] >= max_symb) || (sq->seq[j][i] < 0)) {
        char *str =
          mprintf (NULL, 0, "Wrong symbol \'%d\' in sequence %d at Pos. %d;\
                            Should be within [0..%d]\n",
                   sq->seq[j][i], j + 1, i + 1, max_symb - 1);
        mes_prot (str);
        m_free (str);
        return (-1);
      }
    }
  }
  return 0;
#undef CUR_PROC
}                               /* sequence_check */

/*============================================================================*/
int sequence_best_model (model ** mo, int model_number, int *sequence,
                         int seq_len, double *log_p)
{
# define CUR_PROC "seqence_best_model"
  double log_ptmp;
  int model_index, i;
  *log_p = -DBL_MAX;
  model_index = -1;
  for (i = 0; i < model_number; i++) {
    foba_logp (mo[i], sequence, seq_len, &log_ptmp);
    if (log_ptmp != +1 && log_ptmp > *log_p) {
      *log_p = log_ptmp;
      model_index = i;
    }
  }
  if (*log_p == -DBL_MAX)
    *log_p = +1;
  return (model_index);
# undef CUR_PROC
}                               /* sequence_best_model */


/*============================================================================*/
void sequence_print (FILE * file, sequence_t * sq)
{
  int i, j;
  fprintf (file, "SEQ = {\n\tO = {\n");
  for (i = 0; i < sq->seq_number; i++) {
    if (sq->seq_id[i] != -1.0)
      fprintf (file, "\t(%10.0f)", sq->seq_id[i]);
    if (sq->seq_label[i] != -1)
      fprintf (file, "\t<%ld>", sq->seq_label[i]);
    if (sq->seq_w[i] != 1)
      fprintf (file, "\t|%.0f|", sq->seq_w[i]);
    fprintf (file, "\t");
    if (sq->seq_len[i] > 0) {
      fprintf (file, "%d", sq->seq[i][0]);
      for (j = 1; j < sq->seq_len[i]; j++)
        fprintf (file, ", %d", sq->seq[i][j]);
    }
    fprintf (file, ";\n");
  }
  fprintf (file, "\t};\n};\n\n");
}                               /* sequence_print */

/*============================================================================*/
 /**/ void sequence_print_xml (FILE * file, sequence_t * sq)
{
  int i, j;
  /* coding missing */
  fprintf (file, "<Sequences type=\"int\" >\n");
  fprintf (file, " <DiscretePD>\n");
  for (i = 0; i < sq->seq_number; i++) {
    fprintf (file, "  %.0f <Sequence", sq->seq_w[i]);
    if (sq->seq_id[i] != -1.0)
      fprintf (file, " id=\"seq%f\" ", sq->seq_id[i]);
    fprintf (file, ">");
    if (sq->seq_label[i] != -1)
      fprintf (file, "<Label>%ld</Label>", sq->seq_label[i]);
    if (sq->seq_len[i] > 0) {
      fprintf (file, "<!-- Length: %d -->", sq->seq_len[i]);
      for (j = 0; j < sq->seq_len[i]; j++)
        fprintf (file, " %d", sq->seq[i][j]);
    }
    fprintf (file, "  </Sequence>\n");
  }
  fprintf (file, " </DiscretePD>\n");
  fprintf (file, "</Sequences>\n");
}                               /* sequence_print_xml */

void sequence_d_print_xml (FILE * file, sequence_d_t * sq)
{
  int i, j;
  /* coding missing */
  fprintf (file, "<Sequences type=\"int\" >\n");
  fprintf (file, " <DiscretePD>\n");
  for (i = 0; i < sq->seq_number; i++) {
    fprintf (file, "  %.0f <Sequence", sq->seq_w[i]);
    if (sq->seq_id[i] != -1.0)
      fprintf (file, " id=\"seq%f\" ", sq->seq_id[i]);
    fprintf (file, ">");
    if (sq->seq_label[i] != -1)
      fprintf (file, "<Label>%ld</Label>", sq->seq_label[i]);
    if (sq->seq_len[i] > 0) {
      fprintf (file, "<!-- Length: %d -->", sq->seq_len[i]);
      for (j = 0; j < sq->seq_len[i]; j++)
        fprintf (file, " %f", sq->seq[i][j]);
    }
    fprintf (file, "  </Sequence>\n");
  }
  fprintf (file, " </DiscretePD>\n");
  fprintf (file, "</Sequences>\n");
}                               /* sequence_print_xml */

/*============================================================================*/

void sequence_mathematica_print (FILE * file, sequence_t * sq, char *name)
{
  int i;
  fprintf (file, "%s = {\n", name);
  for (i = 0; i < sq->seq_number - 1; i++)
    vector_i_print (file, sq->seq[i], sq->seq_len[i], "{", ",", "},");
  /* no comma after last seq. */
  vector_i_print (file, sq->seq[sq->seq_number - 1],
                  sq->seq_len[sq->seq_number - 1], "{", ",", "}");
  fprintf (file, "};\n");
}                               /* sequence_d_mathematica_print */

/*============================================================================*/

void sequence_d_gnu_print (FILE * file, sequence_d_t * sqd)
{
  int i, j;
  for (i = 0; i < sqd->seq_number; i++) {
    for (j = 0; j < sqd->seq_len[i]; j++)
      fprintf (file, "%.8f\n", sqd->seq[i][j]);
    fprintf (file, "\n\n");
  }
}

/*============================================================================*/
void sequence_d_print (FILE * file, sequence_d_t * sqd, int discrete)
{
  int i, j;
  fprintf (file, "SEQD = {\n\tO = {\n");
  for (i = 0; i < sqd->seq_number; i++) {
    if (sqd->seq_id[i] != -1.0)
      fprintf (file, "\t(%10.0f)", sqd->seq_id[i]);
    if (sqd->seq_label[i] != -1)
      fprintf (file, "\t<%ld>", sqd->seq_label[i]);
    if (sqd->seq_w[i] != 1)
      fprintf (file, "\t|%.0f|", sqd->seq_w[i]);
    fprintf (file, "\t");
    if (sqd->seq_len[i] > 0) {
      if (discrete)
        fprintf (file, "%3.0f", sqd->seq[i][0]);
      else {
        if (sqd->seq[i][0] > 500)
          fprintf (file, "%8.0f", sqd->seq[i][0]);
        else
          fprintf (file, "%8.2f", sqd->seq[i][0]);
      }
      for (j = 1; j < sqd->seq_len[i]; j++) {
        if (discrete)
          fprintf (file, ", %3.0f", sqd->seq[i][j]);
        else {
          if (sqd->seq[i][j] > 500)
            fprintf (file, ", %8.0f", sqd->seq[i][j]);
          else
            fprintf (file, ", %8.2f", sqd->seq[i][j]);
        }
      }
    }
    fprintf (file, ";\n");
  }
  fprintf (file, "\t};\n};\n\n");
}                               /* sequence_print */

/*============================================================================*/

void sequence_d_mathematica_print (FILE * file, sequence_d_t * sqd,
                                   char *name)
{
  int i;
  fprintf (file, "%s = {\n", name);
  for (i = 0; i < sqd->seq_number - 1; i++)
    vector_d_print (file, sqd->seq[i], sqd->seq_len[i], "{", ",", "},");
  /* no comma after last seq. */
  vector_d_print (file, sqd->seq[sqd->seq_number - 1],
                  sqd->seq_len[sqd->seq_number - 1], "{", ",", "}");
  fprintf (file, "};\n");
}                               /* sequence_d_mathematica_print */

/*============================================================================*/
void sequence_clean (sequence_t * sq)
{
  /* keep data, only delete references */
  m_free (sq->seq);
  m_free (sq->seq_len);
  m_free (sq->seq_label);
  m_free (sq->seq_id);
  m_free (sq->seq_w);
  /* m_free(sq->states);*/
  sq->seq_number = 0;
  sq->total_w = 0.0;
}                               /* sequence_clean */

/*============================================================================*/
void sequence_d_clean (sequence_d_t * sqd)
{
  /* keep data, only delete references */
  m_free (sqd->seq);
  m_free (sqd->seq_len);
  m_free (sqd->seq_label);
  m_free (sqd->seq_id);
  m_free (sqd->seq_w);
  sqd->seq_number = 0;
  sqd->total_w = 0.0;
}                               /* sequence_d_clean */

/*============================================================================*/
int sequence_free (sequence_t ** sq)
{
# define CUR_PROC "sequence_free"
  /*int i,j;*/

  mes_check_ptr (sq, return (-1));
  if (!*sq)
    return (0);

  /*for (i= 0;i<(*sq)->seq_number;i++){
     for (j= 0;j<(*sq)->seq_len[i];j++){
     printf("seq[%d][%d] = %d\n",i,j,(*sq)->seq[i][j]);
     }  
     } */

  /* matrix_i_free also takes care of  (*sq)->seq */
  if (matrix_i_free (&(*sq)->seq, (*sq)->seq_number) == -1) {
    printf ("Error in sequence_free !\n");
  }

  /* XXX The allocation of state must be fixed XXX*/
  /*** Added attribute to the sequence_t
  if (&(*sq)->states) { 
    matrix_i_free(&(*sq)->states, (*sq)->seq_number);
   }
  ***/

  m_free ((*sq)->seq_len);
  m_free ((*sq)->seq_label);
  m_free ((*sq)->seq_id);
  m_free ((*sq)->seq_w);

  if ((*sq)->states) {
    matrix_i_free (&(*sq)->states, (*sq)->seq_number);
    /*m_free((*sq)->states); */
  }


  if ((*sq)->state_labels) {
    matrix_i_free (&(*sq)->state_labels, (*sq)->seq_number);
    m_free ((*sq)->state_labels_len);

    /*m_free((*sq)->states); */
  }

  m_free (*sq);
  return 0;
# undef CUR_PROC
}                               /* sequence_free */


/*============================================================================*/
int sequence_d_free (sequence_d_t ** sqd)
{
# define CUR_PROC "sequence_d_free"
  mes_check_ptr (sqd, return (-1));
  /* sequence_d_print(stdout,*sqd,0);*/

  matrix_d_free (&(*sqd)->seq, (*sqd)->seq_number);
  m_free ((*sqd)->seq_len);
  m_free ((*sqd)->seq_label);
  m_free ((*sqd)->seq_id);
  m_free ((*sqd)->seq_w);
  m_free (*sqd);
  return 0;
# undef CUR_PROC
}                               /* sequence_d_free */

/*============================================================================*/
sequence_d_t *sequence_d_create_from_sq (const sequence_t * sq)
{
# define CUR_PROC "sequence_d_create_from_sq"
  int j, i;
  sequence_d_t *sqd = NULL;     /* target seq. array */
  if (!(sqd = sequence_d_calloc (sq->seq_number))) {
    mes_proc ();
    goto STOP;
  }
  for (j = 0; j < sq->seq_number; j++) {
    ARRAY_CALLOC (sqd->seq[j], sq->seq_len[j]);
    for (i = 0; i < sq->seq_len[j]; i++)
      sqd->seq[j][i] = (double) (sq->seq[j][i]);
    sqd->seq_len[j] = sq->seq_len[j];
    sqd->seq_label[j] = sq->seq_label[j];
    sqd->seq_id[j] = sq->seq_id[j];
    sqd->seq_w[j] = sq->seq_w[j];
  }
  sqd->seq_number = sq->seq_number;
  sqd->total_w = sq->total_w;
  return (sqd);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_d_free (&sqd);
  return NULL;
#undef CUR_PROC
}                               /* sequence_d_create_from_sq */

/*============================================================================*/
sequence_t *sequence_create_from_sqd (const sequence_d_t * sqd)
{
# define CUR_PROC "sequence_create_from_sqd"
  int j, i;
  sequence_t *sq = NULL;        /* target seq. array */
  if (!(sq = sequence_calloc (sqd->seq_number))) {
    mes_proc ();
    goto STOP;
  }
  for (j = 0; j < sqd->seq_number; j++) {
    ARRAY_CALLOC (sq->seq[j], sqd->seq_len[j]);
    for (i = 0; i < sqd->seq_len[j]; i++) {
      sq->seq[j][i] = m_int (fabs (sqd->seq[j][i]));
    }
    sq->seq_len[j] = sqd->seq_len[j];
    sq->seq_label[j] = sqd->seq_label[j];
    sq->seq_id[j] = sqd->seq_id[j];
    sq->seq_w[j] = sqd->seq_w[j];
  }
  sq->seq_number = sqd->seq_number;
  sq->total_w = sqd->total_w;
  return (sq);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_free (&sq);
  return NULL;
#undef CUR_PROC
}                               /* sequence_create_from_sqd */

/*============================================================================*/

int sequence_max_len (const sequence_t * sqd)
{
  int i, max_len = 0;
  for (i = 0; i < sqd->seq_number; i++)
    if (max_len < sqd->seq_len[i])
      max_len = sqd->seq_len[i];
  return max_len;
}                               /* sequence_max_len */


/*============================================================================*/

int sequence_d_max_len (const sequence_d_t * sqd)
{
  int i, max_len = 0;
  for (i = 0; i < sqd->seq_number; i++)
    if (max_len < sqd->seq_len[i])
      max_len = sqd->seq_len[i];
  return max_len;
}                               /* sequence_d_max_len */

/*============================================================================*/

sequence_d_t *sequence_d_mean (const sequence_d_t * sqd)
{
# define CUR_PROC "sequence_d_mean"
  int i, j, max_len;
  sequence_d_t *out_sqd = NULL;

  max_len = sequence_d_max_len (sqd);
  if (!(out_sqd = sequence_d_calloc (1))) {
    mes_proc ();
    goto STOP;
  }
  ARRAY_CALLOC (out_sqd->seq[0], max_len);
  out_sqd->seq_len[0] = max_len;

  for (i = 0; i < sqd->seq_number; i++)
    for (j = 0; j < sqd->seq_len[i]; j++)
      out_sqd->seq[0][j] += sqd->seq[i][j];

  for (j = 0; j < max_len; j++)
    out_sqd->seq[0][j] /= sqd->seq_number;

  return out_sqd;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  sequence_d_free (&out_sqd);
  return NULL;
# undef CUR_PROC
}                               /* sequence_d_mean */

/*============================================================================*/

double **sequence_d_scatter_matrix (const sequence_d_t * sqd, int *dim)
{
# define CUR_PROC "sequence_d_scatter_matrix"
  int *count, k, l, i, j;
  double **W, *mean;

  *dim = sequence_d_max_len (sqd);
  if (!(W = matrix_d_alloc (*dim, *dim))) {
    mes_proc ();
    goto STOP;
  }

  /* Mean over all sequences. Individual counts for each dimension */
  ARRAY_CALLOC (mean, *dim);
  ARRAY_CALLOC (count, *dim);
  for (i = 0; i < sqd->seq_number; i++) {
    for (l = 0; l < sqd->seq_len[i]; l++) {
      mean[l] += sqd->seq[i][l];
      count[l]++;
    }
  }
  for (l = 0; l < *dim; l++)
    mean[l] /= count[l];
  /* scatter matrix (upper triangle) */
  for (j = 0; j < sqd->seq_number; j++) {
    for (k = 0; k < *dim; k++) {
      for (l = k; l < *dim; l++) {
        if (sqd->seq_len[j] > l)
          W[k][l] += (sqd->seq[j][k] - mean[k]) * (sqd->seq[j][l] - mean[l]);
      }
    }
  }
  /* norm with counts, set lower triangle */
  for (k = 0; k < *dim; k++) {
    for (l = *dim - 1; l >= 0; l--) {
      if (l >= k)
        W[k][l] /= (double) count[l];
      else
        W[k][l] = W[l][k];
    }
  }
  return W;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  matrix_d_free (&W, *dim);
  return NULL;
# undef CUR_PROC
}                               /* sequence_d_scatter_matrix */

/*============================================================================*/
/* dummy function at the moment */

int sequence_d_class (const double *O, int index, double *osum)
{
#define CUR_PROC "sequence_d_class"

  return 0;
# undef CUR_PROC
}                               /* sequence_d_class */

/*============================================================================*/

/* divide given field of seqs. randomly into two different fields. Also do 
   allocating. train_ratio determines approx. the fraction of seqs. that go
   into the train_set and test_set resp.
*/

int sequence_d_partition (sequence_d_t * sqd, sequence_d_t * sqd_train,
                          sequence_d_t * sqd_test, double train_ratio)
{
#define CUR_PROC "sequence_d_partition"

  double p;
  sequence_d_t *sqd_dummy = NULL;
  int i;
  long total_seqs, cur_number;

  total_seqs = sqd->seq_number;
  if (total_seqs <= 0) {
    mes_prot ("Error: number of seqs. less or equal zero\n");
    goto STOP;
  }
  /* waste of memory but avoids to many reallocations */
  sqd_dummy = sqd_train;
  for (i = 0; i < 2; i++) {
    ARRAY_CALLOC (sqd_dummy->seq, total_seqs);
    ARRAY_CALLOC (sqd_dummy->seq_len, total_seqs);
    ARRAY_CALLOC (sqd_dummy->seq_label, total_seqs);
    ARRAY_CALLOC (sqd_dummy->seq_id, total_seqs);
    ARRAY_CALLOC (sqd_dummy->seq_w, total_seqs);
    sqd_dummy->seq_number = 0;
    sqd_dummy = sqd_test;
  }

  for (i = 0; i < total_seqs; i++) {
    p = GHMM_RNG_UNIFORM (RNG);
    if (p <= train_ratio)
      sqd_dummy = sqd_train;
    else
      sqd_dummy = sqd_test;
    cur_number = sqd_dummy->seq_number;
    ARRAY_CALLOC (sqd_dummy->seq[cur_number], sqd->seq_len[i]);
    /* copy all entries */
    sequence_d_copy_all (sqd_dummy, cur_number, sqd, i);
    sqd_dummy->seq_number++;
  }

  /* reallocs */
  sqd_dummy = sqd_train;
  for (i = 0; i < 2; i++) {
    ARRAY_REALLOC (sqd_dummy->seq, sqd_dummy->seq_number);
    ARRAY_REALLOC (sqd_dummy->seq_len, sqd_dummy->seq_number);
    ARRAY_REALLOC (sqd_dummy->seq_label, sqd_dummy->seq_number);
    ARRAY_REALLOC (sqd_dummy->seq_id, sqd_dummy->seq_number);
    ARRAY_REALLOC (sqd_dummy->seq_w, sqd_dummy->seq_number);
    sqd_dummy = sqd_test;
  }
  return 0;

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return -1;
#undef CUR_PROC
}

/*============================================================================*/

void sequence_d_copy_all (sequence_d_t * target, long t_num,
                          sequence_d_t * source, long s_num)
{

  sequence_d_copy (target->seq[t_num], source->seq[s_num],
                   source->seq_len[s_num]);
  target->seq_len[t_num] = source->seq_len[s_num];
  target->seq_label[t_num] = source->seq_label[s_num];
  target->seq_id[t_num] = source->seq_id[s_num];
  target->seq_w[t_num] = source->seq_w[s_num];
}


/*============================================================================*/
/* Likelihood function in a mixture model:
   sum_k w^k log( sum_c (alpha_c p(O^k | lambda_c)))
*/

int sequence_d_mix_like (smodel ** smo, int smo_number, sequence_d_t * sqd,
                         double *like)
{
#define CUR_PROC "sequence_d_mix_like"
  int i, k, error_seqs = 0;
  double seq_like = 0.0, log_p;

  *like = 0.0;

  for (i = 0; i < sqd->seq_number; i++) {
    seq_like = 0.0;
    for (k = 0; k < smo_number; k++) {
      if (sfoba_logp (smo[k], sqd->seq[i], sqd->seq_len[i], &log_p) != -1) {
        if (log_p > -100)
          seq_like += exp (log_p) * smo[k]->prior;
      }
    }
    /* no model fits */
    if (seq_like == 0.0) {
      error_seqs++;
      *like += (PENALTY_LOGP * sqd->seq_w[i]);
    }
    else
      *like += (log (seq_like) * sqd->seq_w[i]);
  }

  return error_seqs;

#undef CUR_PROC
}                               /* sequence_d_mix_like */
