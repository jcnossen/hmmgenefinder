/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/linkedlist.c
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
*       This file is version $Revision: 1244 $
*                       from $Date: 2005-08-08 21:42:46 +0200 (Mon, 08 Aug 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/

#include "linkedlist.h"
#include <stdlib.h>
#include "mes.h"
#include <ghmm/internal.h>

void i_list_append(i_list * list, int val){
  i_el * last;
  i_el * el;
  el = init_i_el(val);
  if (list->first == NULL) {
    list->first = el;
    list->last = el;
  }
  else {
    last = list->last;
    last->next = el;
    list->last = el;
  }
  list->length++;
}

void i_list_insert(i_list * list, int val) {
  i_el * first;
  i_el * el;
  el = init_i_el(val);
  if (list->first == NULL) {
    list->first = el;
    list->last = el;
  }
  else {
    first = list->first;
    el->next = first;
    list->first = el;
  }
  list->length++;
}

void i_list_print(i_list * list) {
  i_el * el = list->first;
  printf("LIST : ");
  while(el != NULL) {
    printf("%i, ", el->val);
    el = el->next;
  }
  printf("\n");
}

int * i_list_to_array(i_list * list) {
#define CUR_PROC "i_list_to_array"
  int counter = 0;
  int * array;
  i_el * el;
  ARRAY_CALLOC (array, list->length);
  el = list->first;
  while(el != NULL) {
    array[counter] = el->val;
    el = el->next;
    counter++;
  }
  return array;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  m_free(array);
  return NULL;
#undef CUR_PROC
}

i_list * init_i_list() {
#define CUR_PROC "init_i_list"
  i_list * list;

  ARRAY_CALLOC (list, 1);
  list->first = NULL;
  list->last = NULL;
  list->length = 0;
  return list;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  free_i_list(list);
  return NULL;
#undef CUR_PROC
}

i_el * init_i_el(int val) {
#define CUR_PROC "init_i_el"
  i_el * el;
  ARRAY_CALLOC (el, 1);
  el->next = NULL;
  el->val = val;
  return el;

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  free(el);
  return NULL;
#undef CUR_PROC
}

int free_i_list(i_list * list) {
  i_el * el;
  i_el * next;
  el = list->first;
  while(el != NULL) {
    next = el->next;
    free(el);
    el = next;
  }
  return 0;
}
