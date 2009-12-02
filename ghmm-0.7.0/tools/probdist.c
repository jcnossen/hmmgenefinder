/*******************************************************************************
  author       : Alexander Scliep
  filename     : ghmm/tools/probdist.c
  created      : 1999-3-18
  $Id: probdist.c 1034 2004-10-28 15:25:56Z benrich $


   synopsis:    probdist model1.hmm model2.hmm

   options:     -t <int>   Sequence length 
                -s         symmetric 

   description: computes a probabilistic distance between the two 
                models

                cf. Juang & Rabiner "A Probabilistic Distance Measure
                for Hidden Markov Models"

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




*******************************************************************************/

#ifdef WIN32
#  include "win_config.h"
#endif

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <math.h>

#include "ghmm/mes.h"
#include "ghmm/const.h"
#include "ghmm/foba.h"
#include "ghmm/matrix.h"
#include "ghmm/model.h"
#include "ghmm/sequence.h"
#include "ghmm/viterbi.h"
#include "ghmm/rng.h"

#ifdef CMODEL_INCLUDED
#include "ghmm/cfoba.h"
#include "ghmm/cmodel.h"
#endif
#include "ghmm/smodel.h"


int main(int argc, char* argv[]) {
# define CUR_PROC "main"

  int mo_number, smo_number;
#ifdef CMODEL_INCLUDED
  int cmo_number;
#endif
  int discrete = 0, smodelflag = 0;
  int T, i, j;
  model **mo = NULL;
#ifdef CMODEL_INCLUDED
  cmodel **cmo = NULL;
#endif
  smodel **smo = NULL;
  double d;
  /*  double dangle, ddiff, s1, s2; */

  if (argc != 5) {
    printf("Usage: main <Model File> <discrete-flag> <T> <smodel flag>\n");
    exit(1);
  }

  discrete = atoi(argv[2]);
  T = atoi(argv[3]);
  smodelflag = atoi(argv[4]);

  ghmm_rng_init();
  if (smodelflag) {
    smo = smodel_read(argv[1], &smo_number);
    if (!smo)  {mes_proc(); return -1;}  
    if (smo_number < 2) {
      printf("Need at least two HMMs to compare (read %d)\n", smo_number);
      return -1;
    }
    for (i = 0; i < smo_number - 1; i++)
      for (j = i + 1; j < smo_number; j++) {
	printf("#----- mo[%d], mo[%d] \n", i , j);
	/* syntax prob_dist: (smo1, smo2, total seqlen., symmetric, verbose) */
	d = smodel_prob_distance(smo[i], smo[j], T, 1, 0);
	printf("probdist = %f\n",d);

      }

  }  
  else if (discrete) {

    mo = model_read(argv[1], &mo_number);
    if (!mo) {mes_proc(); return -1;}        

    if (mo_number < 2) {
      printf("Need at least two HMMs to compare\n");
      return -1;
    }

    printf("#----- mo[0], mo[1] \n");
    d = model_prob_distance(mo[0],mo[1], T, 0, 1);
    printf("d=%f\n",d);
    
    printf("#----- mo[1], mo[0] \n");
    d = model_prob_distance(mo[1],mo[0], T, 0, 1);
    printf("d=%f\n",d);
    
    printf("#----- mo[0], mo[1] \n");
    d = model_prob_distance(mo[0],mo[1], T, 0, 0);
    printf("d=%f\n",d);
    
    printf("#----- mo[1], mo[0] \n");
    d = model_prob_distance(mo[1],mo[0], T, 0, 0);
    printf("d=%f\n",d);
    
    printf("#----- mo[0], mo[1] \n");
    d = model_prob_distance(mo[0],mo[1], T, 1, 1);
    printf("d=%f\n",d);
    
    printf("#----- mo[0], mo[1] \n");
    d = model_prob_distance(mo[0],mo[1], T, 1, 0);
    printf("d=%f\n",d);
    
  }    

#ifdef CMODEL_INCLUDED
  else {
    
    cmo = cmodel_read(argv[1], &cmo_number);
    if (!cmo) {mes_proc(); return -1;}        
    
    if (cmo_number < 2) {
      printf("Need at least two CHMMs to compare\n");
      return -1;
    }
    
    printf("#----- cmo[0], cmo[1] \n");
    d = cmodel_prob_distance(cmo[0], cmo[1], T, 0, 1);
    printf("d=%f\n",d);
    
    printf("#----- cmo[1], cmo[0] \n");
    d = cmodel_prob_distance(cmo[1], cmo[0], T, 0, 1);
    printf("d=%f\n",d);
    
    printf("#----- cmo[0], cmo[1] \n");
    d = cmodel_prob_distance(cmo[0], cmo[1], T, 0, 0);
    printf("d=%f\n",d);
    
    printf("#----- cmo[1], cmo[0] \n");
    d = cmodel_prob_distance(cmo[1], cmo[0], T, 0, 0);
    printf("d=%f\n",d);
    
    printf("#----- cmo[0], cmo[1] \n");
    d = cmodel_prob_distance(cmo[0], cmo[1], T, 1, 1);
    printf("d=%f\n",d);
    
    printf("#----- cmo[0], cmo[1] \n");
    d = cmodel_prob_distance(cmo[0], cmo[1], T, 1, 0);
    printf("d=%f\n",d);
    
    /* coemission likelihood */
    printf("#----- cmo[0]/cmo[1] \n");
    if (cmodel_coemission_likelihood(cmo[0], cmo[1], &d) == -1) d = -1;
    printf("Coemission Likelihood = %e\n",d); 
    printf("#----- cmo[0]/cmo[0] \n");
    if (cmodel_coemission_likelihood(cmo[0], cmo[0], &d) == -1) d = -1;
    printf("Coemission Likelihood = %e\n",d); 
    printf("#----- cmo[1]/cmo[1] \n");
    if (cmodel_coemission_likelihood(cmo[1], cmo[1], &d) == -1) d = -1;
    printf("Coemission Likelihood = %e\n",d); 
    
    printf("#----- D_angle, D_diff, S1, S2\n");
    cmodel_measures(cmo[0], cmo[1], &dangle, &ddiff, &s1, &s2);
    printf("D_angle = %e\n", dangle);
    printf("D_diff  = %e\n", ddiff); 
    printf("S1      = %e\n", s1);
    printf("S2      = %e\n", s2); 
  }
#endif /* CMODEL_INCLUDED*/
  
  return 0;
}
