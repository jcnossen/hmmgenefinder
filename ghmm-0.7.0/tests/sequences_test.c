/*******************************************************************************
  author       : Achim Gädke
  filename     : ghmm/tests/sequences_test.c
  created      : DATE: Thu 26. June 2001
  $Id: sequences_test.c 1195 2005-07-05 15:05:37Z schliep $

Copyright (C) 1998-2005 Alexander Schliep
Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
Copyright (C) 2002-2005 Max-Planck-Institut fuer Molekulare Genetik, Berlin

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA




*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ghmm/sequence.h>

void sequence_alloc_print(void)
{
  sequence_t* seq_array;
  int i;

  seq_array= sequence_calloc(1);
  seq_array->seq_len[0]=10;
  seq_array->seq_label[0]=100;
  seq_array->seq_id[0]=101.0;
  seq_array->seq[0]=(int*)malloc(seq_array->seq_len[0]*sizeof(int));

  for (i=0; i<seq_array->seq_len[0]; i++)
    seq_array->seq[0][i]=1;

  sequence_print_xml(stdout,seq_array);

  sequence_free(&seq_array);
}

int main()
{
  sequence_alloc_print();
  return 0;
}
