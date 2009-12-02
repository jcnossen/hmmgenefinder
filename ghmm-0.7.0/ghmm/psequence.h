#ifndef PSEQUENCE_H
#define PSEQUENCE_H
struct psequence {
  /** for each alphabet in model->number_of_alphabets there is one int seq **/
  int ** seq;
  /** number of alphabets (same as in model) **/
  int number_of_alphabets;
  /** for each sequence position there are also double values (e.g) Ka **/
  double ** d_value;
  /** number of continous sequences **/
  int number_of_d_seqs;
  /** length of the sequence **/
  int length;
};

typedef struct psequence psequence;

psequence * init_psequence(int length, int number_of_alphabets, int number_of_d_seqs);

int free_psequence(psequence * seq, int number_of_alphabets, int number_of_d_seqs);

void set_discrete_psequence(psequence * seq_pointer, int index, int * int_seq);

int * get_discrete_psequence(psequence * seq_pointer, int index);

void set_continuous_psequence(psequence * seq_pointer, int index, double * d_seq);

double * get_continuous_psequence(psequence * seq_pointer, int index);

psequence * slice_psequence(psequence * seq_pointer, int start, int stop);

int get_char_psequence(psequence * seq_pointer, int alphabet, int index);

double get_double_psequence(psequence * seq_pointer, int seq_index, int index);
#endif
