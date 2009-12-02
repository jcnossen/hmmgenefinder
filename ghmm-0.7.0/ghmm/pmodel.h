/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: sclass_change.c
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
#ifndef PMODEL_H
#define PMODEL_H
#include "model.h" 
#include "psequence.h"

struct pmodel;

struct pclass_change_context{

    /* Names of class change module/function (for python callback) */
    char* python_module;
    char* python_function;
    
    /** pointer to class function called with seq X, Y and resp indices 
     in the void you can pass the user data */
    int (*get_class)(struct pmodel*, psequence*, psequence*, int, int,void*);
    
    /* space for any data necessary for class switch, USER is RESPONSIBLE */
    void* user_data;
};
typedef struct pclass_change_context pclass_change_context;

struct pstate {
  /** Initial probability */ 
  double pi;
  /** Log of the initial probability */
  double log_pi;
  /** Output probability */
  double *b;
  int order;
  
  /** IDs of the following states */ 
  int *out_id;  
  /** IDs of the previous states */    
  int *in_id;

  /** transition probs to successor states. (for each transition class)*/
  double **out_a; 
  /** transition probs from predecessor states. (for each transition class)*/ 
  double **in_a;
  /** number of transition classes in this state **/
  int kclasses;
  /** pointer to class function   */
  pclass_change_context *class_change; 
  /** int (*get_class)(int*,int); */
  
  /** Transition probability to a successor 
      double *out_a; */
  /** Transition probablity to a precursor 
      double *in_a;*/

  /** Number of successor states */     
  int out_states; 
  /** Number of precursor states */
  int in_states;  
  /** if fix == 1 --> b stays fix during the training */
  int fix;

  int label;
  
  /** EXTENSIONS for Pair HMMs **/
  /** read offset_x characters from sequence X **/
  int offset_x;
  /** read offset_y characters from sequence Y **/
  int offset_y;
  /** which emission alphabet **/
  int alphabet;
};
typedef struct pstate pstate;

/** @name model
    The complete HMM. Contains all parameters, that define a HMM.
*/
struct pmodel {
  /** Number of states */
  int N;
  /** Number of outputs */   
  int M;   
  /** Vector of the states */
  pstate *s; 
  /** The a priori probability for the model.
      A value of -1 indicates that no prior is defined. 
      Note: this is not to be confused with priors on emission
      distributions*/
  double prior;

  /* contains a arbitrary name for the model */
  char* name;
  
   /** Contains bit flags for varios model extensions such as
      kSilentStates, kTiedEmissions (see ghmm.h for a complete list)
  */
  int model_type;

  /** Flag variables for each state indicating whether it is emitting
      or not. 
      Note: silent != NULL iff (model_type & kSilentStates) == 1  */
  int* silent; /*AS*/

  /** Int variable for the maximum level of higher order emissions */
  int maxorder;
  /** saves the history of emissions as int, 
      the nth-last emission is (emission_history * |alphabet|^n+1) % |alphabet|
      see ...*/
  int emission_history;

  /** Flag variables for each state indicating whether the states emissions
      are tied to another state. Groups of tied states are represented
      by their tie group leader (the lowest numbered member of the group).
      
      tied_to[s] == kUntied  : s is not a tied state
      
      tied_to[s] == s        : s is a tie group leader

      tied_to[t] == s        : t is tied to state s

      Note: tied_to != NULL iff (model_type & kTiedEmissions) == 1  */
  int* tied_to; 
  
  /** Note: State store order information of the emissions.
      Classic HMMS have emission order 0, that is the emission probability
      is conditioned only on the state emitting the symbol.

      For higher order emissions, the emission are conditioned on the state s
      as well as the previous emission_order[s] observed symbols.

      The emissions are stored in the state's usual double* b. The order is
      set state.order.

      Note: state.order != NULL iff (model_type & kHigherOrderEmissions) == 1  */
  
  /** background_distributions is a pointer to a
      background_distributions structure, which holds (essentially) an
      array of background distributions (which are just vectors of floating
      point numbers like state.b).

      For each state the array background_id indicates which of the background
      distributions to use in parameter estimation. A value of kNoBackgroundDistribution
      indicates that none should be used.


      Note: background_id != NULL iff (model_type & kHasBackgroundDistributions) == 1  */
  int *background_id;
  background_distributions* bp; 

  /** (WR) added these variables for topological ordering of silent states 
      Condition: topo_order != NULL iff (model_type & kSilentStates) == 1
   */
  int* topo_order; 
  int  topo_order_length;

  /** EXTENSIONS for Pair HMMs **/
  /** total number of alphabets **/
  int number_of_alphabets;
  /** list of sizes of the alphabets **/
  int * size_of_alphabet;
  /** number of double sequences to modify the transition classes */
  int number_of_d_seqs;
  /** maximal offset in sequence X (for the viterbi lookback matrix) **/
  int max_offset_x;
  /** maximal offset in sequence Y (for the viterbi lookback matrix) **/
  int max_offset_y;

  int debug;

};
typedef struct pmodel pmodel;

int pstate_alloc(pstate *state, int M, int in_states, int out_states);

pmodel * init_pmodel();

pclass_change_context * init_pclass_change_context();

void pstate_clean(pstate *my_state);

int pmodel_free(pmodel *mo);

pstate *get_pstateptr(pstate *ary, int index); 

/** functions dealing with the emission indices **/
int pair(int symbol_x, int symbol_y, int alphabet_size, int off_x, int off_y);

int emission_table_size(pmodel * mo, int state_index);
  
int default_transition_class(pmodel * mo, psequence * X, psequence * Y, int index_x, int index_y, void * user_data) ;

void set_to_default_transition_class(pclass_change_context * pccc);

#endif
