/*
  author       : David Posada <dposada@variagenics.com>
  filename     : ghmm/tests/sequences_old_format.c
  created      : DATE: Februar 2002
  $Id: sequences_old_format.c 208 2002-02-24 20:37:26Z achim $
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
  
  
  
*/

/* should be corrected with ghmm/sequence.c version 1.9 */
#include <stdio.h>
#include <ghmm/sequence.h>
   
int main()
{
    int test_result=0;
    const char* double_sequences_file="data/test100.sqd";
    sequence_d_t **sqd = NULL;
    int sqd_number;
    const char* int_sequences_file="data/sequences_old_format.sq";    
    sequence_t **data = NULL;
    int data_number;

    /* read double sequences (this works fine)*/
    fprintf(stderr,"reading double sequences from %s ...",double_sequences_file);
    sqd=sequence_d_read((char*)double_sequences_file, &sqd_number);
    if (sqd==NULL) {
      test_result=1;
      fprintf(stdout, " Failed\n");
    }
    else {
      fprintf(stdout," Done\n");
      sequence_d_free(sqd);
    }


     /* read int sequences (this gives a segmentation fault)*/
    fprintf(stderr,"reading int sequences from %s ...",int_sequences_file);
    data=sequence_read((char*)int_sequences_file,&data_number);
    if (data==NULL) {
      test_result=1;
      fprintf(stdout, " Failed\n");
    }
    else {
      fprintf(stdout," Done\n");
      sequence_free(data);
    }
    return test_result;
}
