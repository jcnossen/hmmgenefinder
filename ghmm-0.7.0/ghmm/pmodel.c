/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/pmodel.c
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
*       This file is version $Revision: 1281 $
*                       from $Date: 2005-08-30 20:36:21 +0200 (Tue, 30 Aug 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/

#include "pmodel.h"
#include "mes.h"
#include <ghmm/internal.h>


int pstate_alloc(pstate * s, int M, int in_states, int out_states) {
# define CUR_PROC "pstate_alloc"
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
  return(res);
# undef CUR_PROC
} /* model_state_alloc */

void pstate_clean(pstate *my_state) {
#define CUR_PROC "pstate_clean"
  if (!my_state) return;
  
  if (my_state->b)
    m_free(my_state->b);
  
  if (my_state->out_id)
    m_free(my_state->out_id);
  
  if (my_state->in_id)
    m_free(my_state->in_id);
  
  if (my_state->out_a)
    m_free(my_state->out_a);
  
  if (my_state->in_a)
    m_free(my_state->in_a);
  
  if (my_state->class_change)
    m_free(my_state->class_change);

  my_state->pi         = 0;
  my_state->b          = NULL;
  my_state->out_id     = NULL;  
  my_state->in_id      = NULL;
  my_state->out_a      = NULL;
  my_state->in_a       = NULL;
  my_state->out_states = 0;
  my_state->in_states  = 0;
  my_state->fix        = 0;  

#undef CUR_PROC
} /* pstate_clean */

/* use this to allocate the memory for a pmodel and set the pointers to NULL */
pmodel * init_pmodel() {
#define CUR_PROC "init_pmodel"
  pmodel * mo;
  ARRAY_CALLOC (mo, 1);

  return mo;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return NULL;
#undef CUR_PROC
}

pclass_change_context * init_pclass_change_context() {
#define CUR_PROC "init_pclass_change_context"
  pclass_change_context * pccc;
  ARRAY_CALLOC (pccc, 1);

  pccc->get_class = &default_transition_class;
  pccc->user_data = NULL;
  return pccc;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return NULL;
#undef CUR_PROC
}

int pmodel_free(pmodel *mo) {
#define CUR_PROC "pmodel_free"
  int i;
  mes_check_ptr(mo, return(-1));
  if( !mo ) return(0);
  if (mo->s) {
    for (i = 0; i < mo->N; i++)
      pstate_clean(&(mo->s[i]));
    m_free(mo->s);
  }
  if (mo->  silent)
    m_free(mo->silent);

  if (mo -> tied_to)
    m_free(mo->tied_to);
  
  /*if ((*mo) -> emission_order) Moved to state
    m_free((*mo)->emission_order);*/
  
  if (mo-> topo_order)
    m_free(mo->topo_order);

  if (mo->number_of_alphabets > 0)
    m_free(mo->size_of_alphabet);
         
  m_free(mo);
  return(0);
#undef CUR_PROC
} /* model_free */  

pstate * get_pstateptr(pstate * ary, int index){ return ary + index; }

/* get the emission index for a pair of symbols 
   if the pair cannot be emmited this returns the size of the emission table 
   plus 1 => emission table[size of table] should be 1 and the actual size
   should be one more... */
int pair(int symbol_x, int symbol_y, int alphabet_size, int off_x, int off_y) {
  if (off_x == 0 && symbol_y >= 0)
    return symbol_y;
  if (off_y == 0 && symbol_x >= 0)
    return symbol_x;
  if (symbol_x < 0 || symbol_y < 0) {
    if (off_x == 0 || off_y == 0) 
      return alphabet_size;
    else
      return (alphabet_size - 1) * (alphabet_size - 1) + 1;
  }
  return(symbol_x * alphabet_size + symbol_y);
}

int emission_table_size(pmodel* mo, int state_index) {
  /* the alphabet is over single sequences so get the maximal index for the
     lookup of emission probabilities and use it to determine the size of
     the lookup table */
  int size =  mo->size_of_alphabet[mo->s[state_index].alphabet];
  return pair(size - 1, size - 1, size, mo->s[state_index].offset_x, mo->s[state_index].offset_y) + 1;
}
 
void print_pstate(pstate * s) {
  int i;

  printf("offset x: %i\n", s->offset_x);
  printf("offset y: %i\n", s->offset_y);
  printf("alphabet: %i\n", s->alphabet);
  printf("kclasses: %i\n", s->kclasses);
  printf("in_state: %i\n", s->in_states);
  for (i=0; i<s->in_states; i++)
    printf("%i ", s->in_id[i]);
  printf("\n");
  printf("probabilities...\n");
}
 
void print_pmodel(pmodel* mo) {
  int i;

  printf("Pair HMM model\n");
  printf("max offset x: %i\n", mo->max_offset_x);
  printf("max offset y: %i\n", mo->max_offset_y);
  printf("Number of states: %i\n", mo->N);
  for (i=0; i<mo->N; i++) {
    printf("State %i:\n", i);
    print_pstate(&(mo->s[i]));
  }
}

int default_transition_class(pmodel * mo, psequence * X, psequence * Y, int index_x, int index_y, void * user_data) {
  return 0;
}

void set_to_default_transition_class(pclass_change_context * pccc) {
  if (pccc){
    pccc->get_class = &default_transition_class;
    pccc->user_data = NULL;
  }
  else
    fprintf(stderr, "set_to_default_transition_class: No class change context\n");
}
