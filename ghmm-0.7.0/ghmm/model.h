/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/model.h
*       Authors:  Benhard Knab, Bernd Wichern, Benjamin Georgi, Alexander Schliep
*
*       Copyright (C) 1998-2004 Alexander Schliep 
*       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
*	Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik, 
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
*       This file is version $Revision: 1356 $ 
*                       from $Date: 2005-09-14 14:01:53 +0200 (Wed, 14 Sep 2005) $
*             last change by $Author: cic99 $.
*
*******************************************************************************/
#ifndef MODEL_H
#define MODEL_H

#ifdef __cplusplus
extern "C" {
#endif

/**@name Discrete HMM model */
/*@{ (Doc++-Group: model) */

/** @name background_distributions
    A container for background distributions to be used in the reestimation. Model
    has an ID (== index) to be used for the arrays background_distributions.order
    and background_distributions.b
*/
  struct background_distributions {
  /** Number of distributions */
    int n;
  /** Number of symbols in alphabet */
    int m;
  /** Order of the respective distribution */
    int *order;
  /** The probabilities */
    double **b;
  };
  typedef struct background_distributions background_distributions;


/** @name state
    The basic structure, keeps all parameters that belong to a state. 
*/
  struct state {
  /** Initial probability */
    double pi;
  /** Output probability */
    double *b;
    int order;

  /** IDs of the following states */
    int *out_id;
  /** IDs of the previous states */
    int *in_id;

  /** transition probabilities to successor states. */
    double *out_a;
  /** transition probabilities from predecessor states. */
    double *in_a;

  /** Number of successor states */
    int out_states;
  /** Number of precursor states */
    int in_states;
  /** if fix == 1 --> b stays fix during the training */
    int fix;

    int label;

  };
  typedef struct state state;


struct coord {
  double x;
  double y;
};
typedef struct coord coord;


/** @name model
    The complete HMM. Contains all parameters, that define a HMM.
*/
  struct model {
  /** Number of states */
    int N;
  /** Number of outputs */
    int M;
  /** Vector of the states */
    state *s;
  /** The a priori probability for the model.
      A value of -1 indicates that no prior is defined. 
      Note: this is not to be confused with priors on emission
      distributions*/
    double prior;

    /* contains a arbitrary name for the model */
    char *name;

   /** Contains bit flags for varios model extensions such as
      kSilentStates, kTiedEmissions (see ghmm.h for a complete list)
  */
    int model_type;

  /** Flag variables for each state indicating whether it is emitting
      or not. 
      Note: silent != NULL iff (model_type & kSilentStates) == 1  */
    int *silent;
      /*AS*/
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
    int *tied_to;

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
    background_distributions *bp;

  /** (WR) added these variables for topological ordering of silent states 
      Condition: topo_order != NULL iff (model_type & kSilentStates) == 1
   */
    int *topo_order;
    int topo_order_length;

  /** pow_lookup is a array of precomputed powers

      It contains in the i-th entry M (alphabet size) to the power of i
      The last entry is maxorder+1
  */
    int *pow_lookup;
    
    
    
    /*  storage of model representation information ( read in from XML ) */
    
    /* emission alphabet  */ 
    char*** alphabet;
    
    /* sizes of the different alphabets (for pair HMMs there might be more than one)  */
    int*  alphabet_size;
    
    /* number of entries in alphabet_size */
    int S;
    
    /* an arry of positions of states for graphical representation */ 
    coord *position;
    
    /* state label alphabet (only for labelled HMMs)  */ 
    char** label_alphabet;
 
    /* size of the label_alphabet */
    int label_size;

  };
  typedef struct model model;
  
  
  


/** @name model_direct
    The complete HMM. Keeps the model parameters in a matrix form. 
*/
  struct model_direct {
  /** Number of states */
    int N;
  /** Number of outputs */
    int M;
  /** Prior for the a priori probability for the model.
      Gets the value -1 if no prior defined. */
    double prior;
  /** Transition matrix  */
    double **A;
  /** Output matrix */
    double **B;
  /** Initial matrix */
    double *Pi;
  /** A vector to know the states where the output should not be trained.
      Default value is 0 for all states. */
    int *fix_state;

    /* XXX additional struct members not addedd here; model_direct is 
       depreciated anyways. Not used by C++ interface */

  };
  typedef struct model_direct model_direct;

/** @name hmm_check_t
    Checks the consistence of the model
  */
  struct hmm_check_t {
  /** Number of rows in the A matrix */
    int r_a;
  /** Number of columns in the A matrix */
    int c_a;
  /** Number of rows in the B matrix */
    int r_b;
  /** Number of columns in the B matrix */
    int c_b;
  /** Length of the phi vector */
    int len_pi;
  /** Length of the fix vector */
    int len_fix;
  };
  typedef struct hmm_check_t hmm_check_t;

#ifdef __cplusplus
}
#endif
/*
  Important: The inclusion of sequence.h ist not done before this point in
  order to avoid error by compiling.
*/
#include <ghmm/sequence.h>
#include <ghmm/scanner.h>
#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------------*/
/** 
    binary algorithm to compute powers of integers efficiently
    see Knuth, TAOCP, Vol 2, 4.6.3 
    uses if appropiate lookup table from struct model */
   int model_ipow (const model * mo, int x, unsigned int n);

/** Frees the memory of a model.
    @return 0 for succes; -1 for error
    @param mo:  pointer to a model */
  int model_free (model ** mo);

/**
   Reads in ASCII data to initialize an array of models. Memory allocation for
   the models is done here.
   @return array of pointers to the models
   @param filename:   the ASCII input file
   @param mo_number:  filled with number of models read */
  model **model_read (char *filename, int *mo_number);

/**
   Reads in a model, where the model parameters are explicit given in
   matrix form. Memory allocation for the model is also done here.
   @return pointer to the model
   @param s:       scanner
   @param multip:  multiplicity; gives how many copies should 
   be made of the model */
  model *model_direct_read (scanner_t * s, int *multip);

/**
   Produces simple left-right models given sequences. 
   The function "model_generate_from_sequence" is called for each 
   model that should be made. The sequences are read in from the
   ASCII file and thrown away again when leaving the function.
   @return vector of models
   @param s:          scanner
   @param new_models: number of models to produce */
  model **model_from_sequence_ascii (scanner_t * s, long *mo_number);

/** 
    Produces simple left-right models given sequences. The sequences
    are not read in from file, but exists already as a structur.
    @return vector of models
    @param s:          scanner
    @param new_models: number of models to produce */
  model **model_from_sequence (sequence_t * sq, long *mo_number);

/**
   Copies a given model. Allocates the necessary memory.
   @return copy of the model
   @param mo:  model to copy */
  model *model_copy (const model * mo);

/**
   Tests if all standardization requirements of model are fulfilled. 
   (That is, if the sum of the probabilities is 1).
   @return 0 for success; -1 for error
   @param mo:  model to test */
  int model_check (const model * mo);

/**
   Tests if number of states and number of outputs in the models match.
   @return 0 for succes; -1 for error
   @param mo:           vector of models
   @param model_number: numbr of models */
  int model_check_compatibility (model ** mo, int model_number);

/**
   Test if to models are compatible. That means their states and outputs match.
   @return 0 for succes; -1 for error
   @param mo:     first model
   @param m2:     second model */
int model_check_compatibel_models (const model * mo, const model * m2);

/**
   Produces a model, which generates the given sequence with probability 1.
   The model is a strict left-right model with one state for each element 
   in the sequence and the output in state i is the i-th value in the sequence 
   with probability 1. The model also has a final state, a state with no output.
   @return         pointer to the produced model 
   @param seq:      sequence
   @param seq_len:  length of the sequence
   @param anz_symb: number of symbols in the sequence
*/
  model *model_generate_from_sequence (const int *seq, int seq_len,
                                       int anz_symb);

/** 
    Produces sequences to a given model. All memory that is needed for the 
    sequences is allocated inside the function. It is possible to define
    the length of the sequences global (global_len > 0) or it can be set 
    inside the function, when a final state in the model is reach (a state
    with no output). If the model has no final state, the sequences will
    have length MAX_SEQ_LEN.
    @return             pointer to an array of sequences
    @param mo:          model
    @param seed:        initial parameter for the random value generator
                        (an integer). If seed == 0, then the random value
			generator is not initialized.
    @param global_len:  length of sequences (=0: automatically via final states)
    @param seq_number:  number of sequences
*/
  sequence_t *model_generate_sequences (model * mo, int seed, int global_len,
                                        long seq_number, int Tmax);

/**
   Calculates the sum log( P( O | lambda ) ).
   Sequences, that can not be generated from the given model, are neglected.
   @return    log(P)
   @param mo model
   @param sq sequences       
*/
  double model_likelihood (model * mo, sequence_t * sq);


/**
    Set transition from state 'i' to state 'j' to value 'prob'.
    NOTE: No internal checks - model might get broken if used carelessly.
    @param mo model
    @param i state index
    @param j state index
    @param prob probabilitys
    
*/
  void model_set_transition (model * mo, int i, int j, double prob);

/**
   Writes a model in matrix format.
   @param file: output file
   @param mo:   model
*/
  void model_print (FILE * file, model * mo);

/**
   Writes transition matrix of a model.
   @param file: output file
   @param mo:   model
   @param tab:  format: leading tabs
   @param separator: format: seperator for columns
   @param ending:    format: end of a row  
*/
  void model_A_print (FILE * file, model * mo, char *tab, char *separator,
                      char *ending);
/**
   Writes output matrix of a model.
   @param file: output file
   @param mo:   model
   @param tab:  format: leading tabs
   @param separator: format: seperator for columns
   @param ending:    format: end of a row  
*/
  void model_B_print (FILE * file, model * mo, char *tab, char *separator,
                      char *ending);
/**
   Writes initial allocation vector of a matrix.
   @param file: output file
   @param mo:   model
   @param tab:  format: leading Tabs
   @param separator: format: seperator for columns
   @param ending:    format: end of a row  
*/
  void model_Pi_print (FILE * file, model * mo, char *tab, char *separator,
                       char *ending);
/**
   Writes fix vector of a matrix.
   @param file: output file
   @param mo:   model
   @param tab:  format: leading Tabs
   @param separator: format: seperator for columns
   @param ending:    format: end of a row  
*/
  void model_fix_print (FILE * file, model * mo, char *tab, char *separator,
                        char *ending);
/**
   Writes transposed transition matrix of a model.
   @param file: output file
   @param mo:   model
   @param tab:  format: leading Tabs
   @param separator: format: seperator for columns
   @param ending:    format: end of a row  
*/
  void model_A_print_transp (FILE * file, model * mo, char *tab,
                             char *separator, char *ending);
/**
   Writes transposed output matrix of a model.
   @param file: output file
   @param mo:   model
   @param tab:  format: leading Tabs
   @param separator: format: seperator for columns
   @param ending:    format: end of a row  
*/
  void model_B_print_transp (FILE * file, model * mo, char *tab,
                             char *separator, char *ending);

/**
   Writes transposed initial allocation vector of a matrix.
   @param file: output file
   @param mo:   model
   @param tab:  format: leading Tabs
   @param separator: format: seperator for columns
   @param ending:    format: end of a row  
*/
  void model_Pi_print_transp (FILE * file, model * mo, char *tab,
                              char *ending);

/**
   Writes a HMM in matrix format. The input model must be of type
   model_direct, that is, have the parameters saved as matrices.
   @param file:   output file
   @param mo_d:   model of type model_direct
   @param multip: number of copies to write
*/
  void model_direct_print (FILE * file, model_direct * mo_d, int multip);

/** 
    Writes the parameters of a model sorted by states. 
    Is not very concise.   
    @param file: output file
    @param mo:   model
*/
  void model_states_print (FILE * file, model * mo);

/** 
    Frees all memory from a model, sets the pointers to NULL and 
    variables to zero.
    @param mo_d  HMM structure (\Ref{struct model_direct})
    @param check Check structure (\Ref{struct hmm_check_t})
*/
  void model_direct_clean (model_direct * mo_d, hmm_check_t * check);

/** 
    Tests compatibility of the model components.
    @return 0 for success; -1 for failure 
    @param mo_d  HMM structure  (\Ref{struct model_direct})
    @param check Check structure  (\Ref{struct hmm_check_t})
*/
  int model_direct_check_data (model_direct * mo_d, hmm_check_t * check);

/** Computes probabilistic distance of two models
    @return the distance
    @param m0  model used for generating random output
    @param m  model to compare with
    @param maxT  maximum output length (for HMMs with absorbing states multiple
                 sequences with a toal langth of at least maxT will be 
		 generated)
    @param symmetric  flag, whether to symmetrize distance (not implemented yet)
    @param verbose  flag, whether to monitor distance in 40 steps.
                    Prints to stdout (yuk!)
*/
  double model_prob_distance (model * m0, model * m, int maxT, int symmetric,
                              int verbose);

/** 
    Frees all memory from a state, sets the pointers to NULL and 
    variables to zero.
    @author Peter Pipenbacher
    @param my_state  state to clean (\Ref{struct state})
*/
  void state_clean (state * my_state);


  sequence_t *model_label_generate_sequences (model * mo, int seed,
                                              int global_len, long seq_number,
                                              int Tmax);


/** 
	Calculates the right index for emission array b of state j in model mo
	given an observation obs and taking the state order into account,
	returns -1 if state order exceeds number of so far emitted characters
	@param  mo:  model
	@param   j:  state id 
	@param obs:  integer observation to be updated with
	@param   t:  position of obs in sequence (time)
*/
   int get_emission_index (model * mo, int j, int obs, int t);

/**
	Updates emission history of model mo, discarding the oldest and 'adding' the
	new observation by using modulo and multiplication	
	@param  mo:  model to be updated
	@param obs:  integer observation to be updated with
*/
   void update_emission_history (model * mo, int obs);

/**
	Updates emission history of model mo for backward algorithm by 'adding'
	observation obs to the left,
	(example: obs = 3
	          2 0 0 1 --> 3 2 0 0 )
	@param  mo:  model to be updated
	@param obs:  integer observation to be updated with
*/
   void update_emission_history_front (model * mo, int obs);


/**
    Uses vector_normalize in vector.h
    Normalizes the transition and output probs for each state
    in the given model
    @return 0 if normalization went through
    @param mo: model to be normalized
*/
  int model_normalize (model * mo);

/**
   Add a specific level of noise to the model parameters
   @return     :        -1 on error, 0 else
   @param mo   :        a pointer to the model
   @param level:        amount of noise to use,
                        a noise level of 0.0 doesn't change the model
   @param seed :        seed for ramdom number generator
*/
  int model_add_noise (model * mo, double level, int seed);


/** 
   Allocates a new background_distributions struct and assigs the arguments to
   the respective fields. Note: The arguments need allocation outside of this
   function.
   
   @return     :               0 on success, -1 on error
   @param mo   :               one model
   @param cur  :               a id of a state
   @param times:               number of times the state cur is at least evaluated
*/
  int model_apply_duration (model * mo, int cur, int times);


/**
   Apply the background distributions to the emission probabilities of states of
   the model which have one specified (background_id[state_id] != kNoBackgroundDistribution).

   @return    :                -1 on error, 0 else
   @param mo  :                a pointer to the model
   @param background_weight:   a parameter controlling the weight given to the
                               background. Note, should be between 0 and 1.
*/
  int model_apply_background (model * mo, double *background_weight);


/** 
   Allocates a new background_distributions struct and assigs the arguments to
   the respective fields. Note: The arguments need allocation outside of this
   function.
   
   @return    :               new pointer to a background_distributions struct
   @param n   :               number of distributions
   @param order:              orders of the distribtions
   @param B:                  matrix of distribution parameters
*/
  background_distributions *model_alloc_background_distributions (int n,
                                                                  int m,
                                                                  int *orders,
                                                                  double **B);

  background_distributions
    *model_copy_background_distributions (background_distributions * bg);
  int model_free_background_distributions (background_distributions * bg);

/**
   Calculates the background distribution for a sequence_t
   the calculated distribution is uniform, every state with the same order
   has the same background distribution.
   Set more sophisticated background distribution from python.
   Caution: overwrites the pointer to previous defined background distributions!
   @return    : -1 on error, 0 else
   @param mo  :	a pointer to the model
   @param sq  : a pointer to a sequence_t struct
   
*/
  int model_get_uniform_background (model * mo, sequence_t * sq);

/**
   Calculates the squared distance between two compatible models.
   The distance is normalized by the number of parameters of the models.
   @return:      normalized squared distance
   @param mo:    first model
   @param m2:    second model
*/
  double model_distance(const model * mo, const model * m2);
/**
   Copies a given state. Allocates the necessary memory.
   @author Peter Pipenbacher
   @return copy of the state
   @param my_state:  state to copy */
#if 0
  state *state_copy (state * my_state);
#endif

/**
   Copies a given state to a given destination.
   @author Peter Pipenbacher
   @param source:  state to copy
   @param dest:    destination */
#if 0
  void state_copy_to (state * source, state * dest);
#endif

#ifdef __cplusplus
}
#endif
#endif
/*@} (Doc++-Group: model) */
