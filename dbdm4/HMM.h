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

	mvec<Edge> edges;
	string name;
	mvec<float> emissions;
};

class HMM {
public:
	HMM() {}
	~HMM();

	// add a state with nucleotide emissions set to 0
	HMMState* AddState(string name) { return AddState(name, mvec<float>()); }
	HMMState* AddState(string name, const mvec<float>& emission);

	void NormalizeProbabilities();

	mvec<HMMState*> states;

	// setup GHMM model structure
	void BuildModel();
	void TestModel();

	static mvec<const char*> GetGenicCodons();
	static mvec<const char*> GetStartCodons();
	static mvec<const char*> GetStopCodons();

	void AddStopCodons();
	void AddStartCodons();

protected:
	GHMM_Model *ghmm_mdl; // ghmm model structure
};

