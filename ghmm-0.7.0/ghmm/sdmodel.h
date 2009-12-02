/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/sdmodel.h
*       Authors:  Wasinee Rungsarityotin, Benjamin Georgi
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
*       This file is version $Revision: 1191 $ 
*                       from $Date: 2005-06-21 11:56:12 +0200 (Tue, 21 Jun 2005) $
*             last change by $Author: cic99 $.
*
*******************************************************************************/
#ifndef SDMODEL_H
#define SDMODEL_H

#ifdef __cplusplus
extern "C" {
#endif

/**@name HMM-Modell */
/*@{ (Doc++-Group: model) */

/** @name state
    The basic structure, keeps all parameters that belong to a state. 
*/
  struct sdstate {
  /** Initial probability */
    double pi;
  /** Output probability */
    double *b;
  /** ID of the following state */
    int *out_id;
  /** ID of the previous state */
    int *in_id;

  /** transition probs to successor states. It is a
   matrix in case of mult. transition matrices (COS > 1)*/
    double **out_a;
  /** transition probs from predecessor states. It is a
   matrix in case of mult. transition matrices (COS > 1) */
    double **in_a;

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
    char *label;
    /* XXX Specific variable for ProfileHMM to count the number of
       match states. Not elegant solution.
       WS: if 1 then counts me, 0 don't count me */
    int countme;
  };
  typedef struct sdstate sdstate;

/** @name model
    The complete HMM. Contains all parameters, that define a HMM.
*/
  struct sdmodel {
  /** Number of states */
    int N;
  /** Number of outputs */
    int M;
 /** smodel includes continuous model with one transition matrix 
      (cos  is set to 1) and an extension for models with several matrices
      (cos is set to a positive integer value > 1).*/
    int cos;
  /** Vector of the states */
    sdstate *s;
  /** Prior for the a priori probability for the model.
      A value of -1 indicates that no prior is defined. */
    double prior;

  /** Contains bit flags for various model extensions such as
      kSilentStates, kTiedEmissions (see ghmm.h for a complete list)
  */

  /** pointer to class function   */
    int (*get_class) (int *, int);

    /*int (*get_class)(const double*,int,double*); */

  /** Contains bit flags for various model extensions such as
      kSilentStates, kTiedEmissions (see ghmm.h for a complete list)
  */
    int model_type;

  /** Flag variables for each state indicating whether it is emitting
      or not. 
      Note: silent != NULL iff (model_type & kSilentStates) == 1  */
    int *silent;
      /*AS*/ int topo_order_length;
      /*WR*/ int *topo_order;
    /*WR*/};
  typedef struct sdmodel sdmodel;


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


  /** Frees the memory of a model.
      @return 0 for succes; -1 for error
      @param mo:  pointer to a model */
  int sdmodel_free (sdmodel ** mo);

  int sdmodel_initSilentStates (sdmodel * mo);

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
	  @param T_max:  maximal number of consecutive silent states in model (used to
	  identify silent circles).
  */
  sequence_t *sdmodel_generate_sequences (sdmodel * mo, int seed,
                                          int global_len, long seq_number,
                                          int Tmax);


  /**
     Copies a given model. Allocates the necessary memory.
     @return copy of the model
     @param mo:  model to copy */
  sdmodel *sdmodel_copy (const sdmodel * mo);

  /** Utility for converting between single discrete model and switching model */
  model *sdmodel_to_model (const sdmodel * mo, int kclass);

  /** */
  void model_to_sdmodel (const model * mo, sdmodel * smo, int klass);

  /**
     Writes a model in matrix format.
     @param file: output file
     @param mo:   model
  */
  void sdmodel_print (FILE * file, sdmodel * mo);


  /**
     Writes transition matrix of a model.
     @param file: output file
     @param mo:   model
     @param tab:  format: leading tabs
     @param separator: format: seperator for columns
     @param ending:    format: end of a row  
  */
  void sdmodel_Ak_print (FILE * file, sdmodel * mo, int k, char *tab,
                         char *separator, char *ending);
  /**
     Writes output matrix of a model.
     @param file: output file
     @param mo:   model
     @param tab:  format: leading tabs
     @param separator: format: seperator for columns
     @param ending:    format: end of a row  
  */
  void sdmodel_B_print (FILE * file, sdmodel * mo, char *tab, char *separator,
                        char *ending);

  /**
     Writes initial allocation vector of a matrix.
     @param file: output file
     @param mo:   model
     @param tab:  format: leading Tabs
     @param separator: format: seperator for columns
     @param ending:    format: end of a row  
  */
  void sdmodel_Pi_print (FILE * file, sdmodel * mo, char *tab,
                         char *separator, char *ending);

  /*============================================================================*/
  /** sdviterbi is working for switching discrete model
   *  sdmodel_topo_ordering -- need to be implemented with DFS (as in model_util.c)
   *============================================================================
   **/

  void sdmodel_topo_ordering (sdmodel * mo);

  int *sdviterbi (sdmodel * mo, int *o, int len, double *log_p);

  /** Forward-Algorithm.
      Calculates alpha[t][i], scaling factors scale[t] and log( P(O|lambda) ) for
      a given double sequence and a given model.
      @param smo      model
      @param O        sequence
      @param length: length of sequence
      @param alpha:  alpha[t][i]
      @param scale:   a reference for double type, scale factors
      @param log\_p:  a reference for double type, log likelihood log( P(O|lambda) )
      @return 0 for success, -1 for error
  */
  int sdfoba_forward (sdmodel * mo, const int *O, int len, double **alpha,
                      double *scale, double *log_p);


  /** Descale
      descales the alpha matrix from the forward algorithm
      @param alpha: alpha matrix from forward
      @param scale: scale vector from forward
      @param t:     number of timesteps
      @param n:     number of states
      @param newalpha: unscaled alpha matrix
      @return 0 for success, -1 for error
  */
  int sdfoba_descale (double **alpha, double *scale, int t, int n,
                      double **newalpha);


/**
   Calculates the sum log( P( O | lambda ) ).
   Sequences, that can not be generated from the given model, are neglected.
   @return    log(P)
   @param mo model
   @param sq sequences       
*/
  double sdmodel_likelihood (sdmodel * mo, sequence_t * sq);


/** 
    Writes the parameters of a model sorted by states. 
    Is not very concise.   
    @param file: output file
    @param mo:   model
*/
  void sdmodel_states_print (FILE * file, sdmodel * mo);


#ifdef __cplusplus
}
#endif
#endif
/*@} (Doc++-Group: model) */
