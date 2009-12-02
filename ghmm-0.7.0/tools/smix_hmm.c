/*******************************************************************************
  author       : Bernd Wichern
  filename     : /zpr/bspk/src/hmm/ghmm/tools/smix_hmm.c
  created      : TIME: 17:23:38     DATE: Tue 18. September 2001
  last-modified: TIME: 13:46:48     DATE: Wed 19. September 2001
  $Id: smix_hmm.c 1034 2004-10-28 15:25:56Z benrich $
*******************************************************************************/

#ifdef WIN32
#  include "win_config.h"
#endif

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ghmm/mes.h>
#include <ghmm/rng.h>
#include <ghmm/smodel.h>
#include <ghmm/sequence.h>
#include <ghmm/smixturehmm.h>
#include <ghmm/matrix.h>

static int smix_hmm_run(int argc, char* argv[]) {
#define CUR_PROC "smix_hmm_run"
  int k, exitcode = -1, smo_number, sqd_fields;
  sequence_d_t **sqd = NULL;
  smodel **smo = NULL;
  double **cp = NULL;
  FILE *outfile = NULL;

  /* read sequences and initial models */
  sqd = sequence_d_read(argv[1], &sqd_fields);
  if (!sqd) {mes_proc(); goto STOP;}
  if (sqd_fields > 1)
    printf("Warning: Seq. File contains multiple Seq. Fields; use only the first one\n");
  smo = smodel_read(argv[2], &smo_number);
  if (!smo) {mes_proc(); goto STOP;}

  /* open output file */
  if(!(outfile = mes_fopen(argv[3], "wt"))) {mes_proc(); goto STOP;}
  
  /* matrix for component probs., */
  cp = matrix_d_alloc(sqd[0]->seq_number, smo_number);
  if (!cp) { mes_proc(); goto STOP;}

  /* set last arg in smixturehmm_init() : 
     1 = strict random partition; cp = 0/1
     2. smap_bayes from initial models
     3. cp = 1 for best model, cp = 0 for other models 
     4. open
     5. no start partition == equal cp for each model
  */
  if (smixturehmm_init(cp, sqd[0], smo, smo_number, 5) == -1) {
    mes_proc(); goto STOP;
  }
  /* clustering */
  if (smixturehmm_cluster(outfile, cp, sqd[0], smo, smo_number) == -1) {
    mes_proc(); goto STOP;
  }

  /* print trained models */
  for (k = 0; k < smo_number; k++)
    smodel_print(outfile, smo[k]);  

  if (outfile) fclose(outfile);
  exitcode = 0;
 STOP:
  return exitcode;
# undef CUR_PROC
}


/*============================================================================*/
int main(int argc, char* argv[]) {
#define CUR_PROC "smix_hmm_main"
  int exitcode = -1;

  if (argc != 4 && argc != 5) {
    printf("Insufficient arguments. Usage: \n");
    printf("mix_hmm [Seq.File] [InitModel File] [Out File] <seed>\n");
    goto STOP;
  }
  ghmm_rng_init();
  if (argc == 5)
      GHMM_RNG_SET(RNG,atoi(argv[4]));
    else {      
      ghmm_rng_timeseed(RNG); 
    }
  
  exitcode = smix_hmm_run(argc, argv);
  /*------------------------------------------------------------------------*/
 STOP:
  mes(MES_WIN, "\n(%2.2T): Program finished with exitcode %d.\n", exitcode );
  mes_exit();
  return(exitcode);
# undef CUR_PROC
} /* main */

   
