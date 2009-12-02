/*******************************************************************************
  author       : Bernd Wichern
  filename     : /zpr/bspk/src/hmm/ghmm/tests/shmm_viterbi_test.c
  created      : TIME: 12:52:37     DATE: Tue 05. June 2001
  last-modified: TIME: 15:33:06     DATE: Tue 05. June 2001
  $Id: shmm_viterbi_test.c 201 2002-02-22 23:48:44Z achim $
*******************************************************************************/

#include <string.h>
#include <ghmm/smodel.h>
#include <ghmm/sequence.h>
#include <ghmm/sviterbi.h>
#include <ghmm/mes.h>

static int viterbi_test(char* argv[]);

/*============================================================================*/

static int viterbi_test(char* argv[]) {
#define CUR_PROC "viterbi_test"
  FILE *outfile = NULL;
  sequence_d_t **sqd = NULL;
  smodel** smo = NULL;
  char outfilename[256];
  int smo_number, sqd_number, *state_seq = NULL, model, t, cnt = 0;
  long i, j;
  double log_p;


  /* read array of models */
  smo = smodel_read(argv[1], &smo_number);
  if (!(smo)) {mes_proc(); goto STOP;}
  
  /* read sequences */
  sqd = sequence_d_read(argv[2], &sqd_number);
  if (!sqd) {mes_proc(); goto STOP;}
  
  strcpy(outfilename, argv[3]);
  if(!(outfile = mes_fopen(outfilename, "wt"))) {mes_proc(); goto STOP;}

  /* calculate viterbi path for every possible sequence-model combination */
  for (model = 0; model < smo_number; model++) {
    for (i = 0;  i < sqd_number; i++) {
      for (j = 0; j < sqd[i]->seq_number; j++) {
	cnt++;
	state_seq = sviterbi(smo[model], sqd[i]->seq[j], sqd[i]->seq_len[j],
			     &log_p);
	if (state_seq == NULL) {mes_proc(); goto STOP;}
	fprintf(outfile, "%d %ld (Seq.ID %d): logp %.4f\t", model, i, 
		(int) sqd[i]->seq_id[j], log_p);
	for (t = 0; t < sqd[i]->seq_len[j]; t++)
	  fprintf(outfile, "%2d ", state_seq[t]);
	fprintf(outfile, "\n");
	m_free(state_seq);
	if (!cnt%100) printf("%d\n", cnt);
      }
    }
  }


STOP:
  return -1;
# undef CUR_PROC
}

/*============================================================================*/

int main(int argc, char* argv[]) {

  if (argc != 4) {
    printf("Usage: shmm_viterbi_test <Model File> <Sequence File> <Output File> \n");
    exit(1);
  }

  return viterbi_test(argv);
}

