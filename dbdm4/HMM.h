#pragma once

#include "DNAUtil.h"

using std::string;

typedef state GHMM_State;
typedef model GHMM_Model;

class HMMState {
public:
	HMMState(const string& name) : name(name) {}

	struct Edge {
		Edge(float P, HMMState* dst) : prob(P), dst(dst) {}

		float prob;
		HMMState* dst;
	};

	void AddEdge(HMMState* s, float prob);

	mvec<Edge> edges;
	string name;
	mvec<float> emissions;
};

class HMM {
public:
	HMM() { ghmm_mdl = 0; }
	~HMM();

	// add a state with nucleotide emissions set to 0
	HMMState* AddState(string name) { return AddState(name, mvec<float>()); }
	HMMState* AddState(string name, const mvec<float>& emission);

	mvec<HMMState*> states;

	void TestModel();
	// generate a random sequence based on model
	// BuildModel needs to be done first
	mvec<int> GenerateSequence(int len);

	void AddStopCodons();
	void AddStartCodons();

	void MergeHMM(HMM* hmm);

	// setup GHMM model structure
	void BuildModel();

protected:
	void NormalizeProbabilities();

	GHMM_Model *ghmm_mdl;
};

class HMM_Intergenic : public HMM
{

};