/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/psequence.c
*       Authors:  Matthias Heinig
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
*       This file is version $Revision: 1267 $
*                       from $Date: 2005-08-17 19:29:10 +0200 (Wed, 17 Aug 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/

#include "psequence.h"
#include <stdlib.h>
#include "mes.h"
#include "matrix.h"
#include <ghmm/internal.h>

psequence * init_psequence(int length, int number_of_alphabets, int number_of_d_seqs) {
#define CUR_PROC "init_psequence"
  psequence * seq;

  ARRAY_MALLOC (seq, 1);

  seq->length = length;
  seq->number_of_alphabets = number_of_alphabets;
  seq->number_of_d_seqs = number_of_d_seqs;
  seq->seq = NULL;
  seq->d_value = NULL;
  if (number_of_alphabets > 0) {
    seq->seq = matrix_i_alloc(number_of_alphabets, length); 
    if (!(seq->seq)) goto STOP;
  }
  if (number_of_d_seqs > 0) {
    seq->d_value = matrix_d_alloc(number_of_d_seqs, length);
    if (!(seq->d_value)) goto STOP;
  }

  return seq;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  free_psequence(seq, number_of_alphabets, number_of_d_seqs);
  return NULL;
#undef CUR_PROC
}

int free_psequence(psequence * seq, int number_of_alphabets, int number_of_d_seqs) {
#define CUR_PROC "free_psequence"
  int i;
  mes_check_ptr(seq, return(-1));
  if ( seq == NULL ) return(0);
  if (seq->seq != NULL) {
    for (i=0; i<number_of_alphabets; i++) 
      m_free(seq->seq[i]);
    m_free(seq->seq);
  }
  if (seq->d_value != NULL) {
    for (i=0; i<number_of_d_seqs; i++)
      m_free(seq->d_value[i]);
    m_free(seq->d_value);
  }
  m_free(seq);
  return 0;
#undef CUR_PROC
}

void set_discrete_psequence(psequence * seq_pointer, int index, int * int_seq) {
  seq_pointer->seq[index] = int_seq;
}
  
void set_continuous_psequence(psequence * seq_pointer, int index, double * d_seq) {
  seq_pointer->d_value[index] = d_seq;
}

int * get_discrete_psequence(psequence * seq_pointer, int index){
  return seq_pointer->seq[index];
}

double * get_continuous_psequence(psequence * seq_pointer, int index){
  return seq_pointer->d_value[index];
}

psequence * slice_psequence(psequence * seq_pointer, int start, int stop){
  int i, j;
  psequence * slice;

  if (stop > seq_pointer->length) {
    fprintf(stderr, "Slice: sequence index (%i) out of bounds (%i)\n", 
	    stop, seq_pointer->length);
  }
  slice = init_psequence(stop - start, seq_pointer->number_of_alphabets,
			  seq_pointer->number_of_d_seqs);
  for (i=start; i<stop; i++){
    for (j=0; j<slice->number_of_alphabets; j++)
      slice->seq[j][i-start] = seq_pointer->seq[j][i];
    for (j=0; j<slice->number_of_d_seqs; j++)
      slice->d_value[j][i-start] = seq_pointer->d_value[j][i];
  }
  return slice;
}

int get_char_psequence(psequence * seq_pointer, int alphabet, int index){
  if (alphabet < seq_pointer->number_of_alphabets) {
    if (index < 0)
      return -1;
    if (index < seq_pointer->length) {
      return seq_pointer->seq[alphabet][index];
    }
    else {
      fprintf(stderr, "index (%i) larger than length (%i)!", index, seq_pointer->length);
      return -1;
    }
  }
  else {
    fprintf(stderr, "alphabet (%i) larger than number of alphabets (%i)!",
	    alphabet, seq_pointer->number_of_alphabets);
    return -1;
  }
}

double get_double_psequence(psequence * seq_pointer, int seq_index, int index){
  if (seq_index < seq_pointer->number_of_d_seqs) {
    if (index < 0)
      return 0;
    if (index < seq_pointer->length) {
      return seq_pointer->d_value[seq_index][index];
    }
    else {
      fprintf(stderr, "index (%i) larger than length (%i)!", index, seq_pointer->length);
      return 0;
    }
  }
  else {
    fprintf(stderr, "seq_index (%i) larger than number of seq_indexs (%i)!",
	    seq_index, seq_pointer->number_of_d_seqs);
    return 0;
  }
} 
