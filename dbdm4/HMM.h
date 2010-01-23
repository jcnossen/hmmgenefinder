/*
HMM class, wraps the GHMM into some easier usable interface
*/
#pragma once

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

	mvec<Edge> inputs;
	mvec<Edge> outputs;
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

	void MergeHMM(HMM* hmm);
	void ListStates();

	// setup GHMM model structure
	void BuildModel();

	// train, and feed back into HMMState's
	void BaumWelch(const mvec<int>& seq);
	mvec<int> ViterbiPath(const mvec<int>& nt);

protected:
	void NormalizeProbabilities();
	void CopyParametersFromModel(); // copy from ghmm_mdl

	GHMM_Model *ghmm_mdl;
};
