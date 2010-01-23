
#include "stdafx.h"

#include "HMM.h"
#include "DNAUtil.h"
#include "Genome.h"

void HMMState::AddEdge( HMMState* s, float prob )
{
	edges.push_back( Edge(prob, s) );
}


HMM::~HMM()
{
	DeleteAll(states);

	if (ghmm_mdl) {
		model_free(&ghmm_mdl);
	}
}

HMMState* HMM::AddState( std::string name, const mvec<float>& emissions )
{
	HMMState* s = new HMMState(name);

	if (!emissions.empty())
		s->emissions = emissions;

	states.push_back(s);
	return s;
}

template<typename T> T* alloc(int n=1) {
	void *d = mes_malloc(sizeof(T) * n);
	memset(d, 0, sizeof(T) * n);
	return (T*)d;
}

// an edge that uses a index instead of HMMState*
struct IndexEdge {
	int index;
	float prob;
};


// Setup GHMM model structure
void HMM::BuildModel()
{
	NormalizeProbabilities();
	
	if (ghmm_mdl) 
		model_free(&ghmm_mdl);

	GHMM_Model * mdl = alloc<GHMM_Model>();
	ghmm_mdl = mdl;

	mdl->M = 4;
	mdl->N = states.size();

	mdl->s = alloc<GHMM_State>(mdl->N);
	mdl->prior = -1;
	mdl->model_type = kSilentStates; // enable silent states
	mdl->silent = alloc<int>(states.size());

	// map states to indices
	std::map< HMMState*, int > stateToIndex;
	for(int i=0;i<states.size();i++)
		stateToIndex[states[i]]=i;

	// generate input lists for each state (only outputs are stored in the edges)
	mvec< mvec< IndexEdge > > stateInputs;
	stateInputs.resize(states.size());
	for(int i=0;i<states.size();i++) {
		for (int j=0;j<states[i]->edges.size();j++) {
			HMMState::Edge& edge = states[i]->edges[j];
			int stateIndex = stateToIndex[edge.dst];
			IndexEdge ie;
			ie.index = i;
			ie.prob = edge.prob;
			stateInputs[stateIndex].push_back(ie);
		}
	}

	// fill state info
	for(int i=0;i<states.size();i++) {
		HMMState* src = states[i];
		GHMM_State* dst = &mdl->s[i];

		dst->b = alloc<double>(mdl->M);
		if (src->emissions.empty()) {
			mdl->silent[i] = 1;
		} else {
			mdl->silent[i] = 0;

			// set state emissions (output)
			for (int j=0;j<mdl->M;j++) 
				dst->b[j] = src->emissions[j];
		}
		dst->pi = 1.0f; // initial probability (TODO) ?

		// Build state inputs
		mvec<IndexEdge> &inputs = stateInputs[i];
		dst->in_states = inputs.size();
		dst->in_a = alloc<double>(inputs.size());
		dst->in_id = alloc<int>(inputs.size());

		// normalize inputs
		float invTotalProb = 1.0f / sum(members(inputs, &IndexEdge::prob));

		for (int j=0;j<inputs.size();j++) {
			dst->in_a[j]=inputs[j].prob * invTotalProb;
			dst->in_id[j]=inputs[j].index;
		}

		// Build state outputs
		dst->out_states = src->edges.size();
		dst->out_a = alloc<double>(dst->out_states);
		dst->out_id = alloc<int>(dst->out_states);
		for (int j =0 ;j<src->edges.size();j++) {
			HMMState::Edge& e = src->edges[j];
			dst->out_a[j] = e.prob;
			dst->out_id[j] = stateToIndex[e.dst];
		}
	}
}

// Normalize all state output probabilities
void HMM::NormalizeProbabilities()
{
	for (int i=0;i<states.size();i++) {
		HMMState& st = *states[i];

		if (!st.emissions.empty())
			st.emissions /= sum( st.emissions );

		float invEdgeProbSum = 1.0f / sum( members(st.edges, &HMMState::Edge::prob) );
		for (int j=0;j<st.edges.size();j++)
			st.edges[j].prob *= invEdgeProbSum;
	}
}

	// 
	// 	% Adds stop codons to the model
	// 		function [] = add_stop(obj)
	// 		id1(1) = obj.add_state('stop_TAA_TGA_1', HMM.get_emission_prob_for_nucleotide('T'));
	// 	count_TAA = obj.stop_codon_stats(nt2int('T'), nt2int('A'), nt2int('A'));
	// 	count_TGA = obj.stop_codon_stats(nt2int('T'), nt2int('G'), nt2int('A'));
	// 	count_both = count_TAA + count_TGA;
	// 	id1(2) = obj.add_state('stop_TAA_TGA_2', (HMM.get_emission_prob_for_nucleotide('A') * count_TAA + HMM.get_emission_prob_for_nucleotide('G') * count_TGA) ./ count_both);
	// 	id1(3) = obj.add_state('stop_TAA_TGA_3', HMM.get_emission_prob_for_nucleotide('A'));
	// 	obj.add_edge(id1(1), id1(2), 1);
	// 	obj.add_edge(id1(2), id1(3), 1);
	// 
	// 	id2(1) = obj.add_state('stop_TAG_1', HMM.get_emission_prob_for_nucleotide('T'));
	// 	id2(2) = obj.add_state('stop_TAG_2', HMM.get_emission_prob_for_nucleotide('A'));
	// 	id2(3) = obj.add_state('stop_TAG_3', HMM.get_emission_prob_for_nucleotide('G'));
	// 
	// 	obj.add_edge(id2(1), id2(2), 1);
	// 	obj.add_edge(id2(2), id2(3), 1);

void HMM::TestModel()
{
	assert(ghmm_mdl);
	fprintf(stdout,"transition matrix:\n");
	model_A_print(stdout,ghmm_mdl,""," ","\n");
	fprintf(stdout,"observation symbol matrix:\n");
	model_B_print(stdout,ghmm_mdl,""," ","\n");

	sequence_t *seq = model_generate_sequences(ghmm_mdl,0,20,1,100);
	sequence_print(stdout, seq);
	sequence_free(&seq);
}

mvec<int> HMM::GenerateSequence( int len )
{
	sequence_t *seq = model_generate_sequences(ghmm_mdl,0,len,1,100);
	mvec<int> r(seq->seq[0], len);
	sequence_free(&seq);
	return r;
}

void HMM::MergeHMM( HMM* hmm )
{
	states &= hmm->states;
	hmm->states.clear();
}

void HMM::ListStates()
{
	d_trace("%d states:\n");
	for (int i=0;i<states.size();i++) {
		d_trace("\tstate [%d]: %s\n", i, states[i]->name.c_str());
	}
}

mvec<int> HMM::ViterbiPath( const mvec<int>& nt )
{
	int *vpath;
	double log_p;

	assert(ghmm_mdl);

	vpath = viterbi(ghmm_mdl,(int*) &nt[0], nt.size(), &log_p);
	// find end of path 
	int pos = 0;
	while (vpath[pos] >= 0)
		pos++;

	mvec<int> path(vpath, vpath + pos);
	m_free(vpath);

	return path;
}

void HMM::Forward()
{
}

