
#include "stdafx.h"

#include "HMM.h"
#include "DNAUtil.h"
#include "Genome.h"
#include "CfgParser.h"

void HMMState::AddEdge( HMMState* s, float prob )
{
	outputs.push_back(Edge(prob, s));
	s->inputs.push_back(Edge(1.0f, this));
}


HMM::~HMM()
{
	DeleteAll(states);

	if (ghmm_mdl) {
		model_free(&ghmm_mdl);
	}
}

HMMState* HMM::AddState( std::string name, const mvec<double>& emissions )
{
	HMMState* s = new HMMState(name);

	if (!emissions.empty())
		s->emissions = emissions;

	states.push_back(s);
	return s;
}

// For allocating GHMM structure memory
template<typename T> T* alloc(int n=1) 
{
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
		dst->pi = initial_state==src ? 1.0f : 0.0f; // initial probability

		// Build state inputs
		mvec<HMMState::Edge> &inputs = src->inputs;
		dst->in_states = inputs.size();
		dst->in_a = alloc<double>(inputs.size());
		dst->in_id = alloc<int>(inputs.size());

		for (int j=0;j<inputs.size();j++) {
			dst->in_a[j]= inputs[j].prob;
			dst->in_id[j]= stateToIndex[ inputs[j].dst ];
		}

		// Build state outputs
		dst->out_states = src->outputs.size();
		dst->out_a = alloc<double>(dst->out_states);
		dst->out_id = alloc<int>(dst->out_states);
		for (int j =0 ;j<src->outputs.size();j++) {
			HMMState::Edge& e = src->outputs[j];
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

		double invEdgeProbSum = 1.0f / sum( members(st.outputs, &HMMState::Edge::prob) );
		for (int j=0;j<st.outputs.size();j++)
			st.outputs[j].prob *= invEdgeProbSum;

		if (st.inputs.empty()) {
			d_trace("warning: state %s has no inputs\n", st.name.c_str());
		} else {
			invEdgeProbSum = 1.0f / sum( members(st.inputs, &HMMState::Edge::prob) );
			for (int j=0;j<st.inputs.size();j++)
				st.inputs[j].prob *= invEdgeProbSum;
		}
	}
}


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


void HMM::BaumWelch( const mvec< mvec<int>* >& sequences )
{
	sequence_t* st = sequence_calloc(sequences.size());

	for (int i=0;i<sequences.size();i++) {
		mvec<int>* seq = sequences[i];
		st->seq[i] = alloc<int>(seq->size());
		std::copy(seq->begin(), seq->end(), st->seq[i]);
		st->seq_len[i] = seq->size();
	}

	int maxStep = 500;
	reestimate_baum_welch_nstep(ghmm_mdl, st, maxStep, EPS_ITER_BW);
	CopyParametersFromModel();

	sequence_free(&st);
}

// copy the GHMM data back into ours
void HMM::CopyParametersFromModel()
{
	for(int i=0;i<ghmm_mdl->N;i++) {
		HMMState *dst = states[i];
		GHMM_State *src = &ghmm_mdl->s[i];

		for (int j=0;j<dst->outputs.size();j++) {
	//		d_trace("%s.outputs[%d] old: %f, new: %f\n", dst->name.c_str(), j, dst->outputs[j].prob, src->out_a[j]);
			dst->outputs[j].prob = src->out_a[j];
		}

		for (int j=0;j<dst->inputs.size();j++) {
//			d_trace("%s.inputs[%d] old: %f, new: %f\n", dst->name.c_str(), j, dst->inputs[j].prob, src->in_a[j]);
			dst->inputs[j].prob = src->in_a[j];
		}
	}
}

void HMM::ParseConfig(std::string file)
{
	CfgList* cfg=CfgList::LoadFile(file.c_str());

//	CfgList* symbols = cfg->GetList("symbols");

	CfgList* stateList = cfg->GetList("states");

	states.resize(stateList->childs.size());
	for (int i=0;i<stateList->childs.size();i++) {
		CfgList* st = (CfgList*) stateList->childs[i]->value;
		states[i] = new HMMState(st->GetLiteral("name", ""));
	}

	int nEdge=0;

	for (int i=0;i<stateList->childs.size();i++) {
		CfgList* st = (CfgList*) stateList->childs[i]->value;

		CfgList* emit = st->GetList("emit");
		bool allzero=true;
		for (int j=0;j<emit->childs.size();j++) {
			states[i]->emissions.push_back( ((CfgNumeric*) emit->childs[j]->value)->value );
			if ( fabs(states[i]->emissions.back()) < 0.0001)
				allzero=false;

		}
		if (allzero)
			states[i]->emissions.clear();

		CfgList* outputs = st->GetList("outputs");
		for (int j=0;j<outputs->childs.size();j++) {
			double prob;
			int index;

			CfgList* outputcfg =  (CfgList*)outputs->childs[j]->value;
			// assume indices for speed sake
			index = (int)((CfgNumeric*)outputcfg->childs[0]->value)->value;
			prob = ((CfgNumeric*)outputcfg->childs[0]->value)->value;

			index --; // matlab eh..
			assert(index >= 0 && index < states.size());
			states[i]->AddEdge(states[index], prob);
			nEdge ++;
		}

	}

	d_trace("Config %s contains %d states and %d edges\n", file.c_str(), states.size(), nEdge);
}

HMMState* HMM::FindState( string name )
{
	for (int i=0;i<states.size();i++)
		if (states[i]->name == name )
			return states[i];

	return 0;
}