
#include "stdafx.h"

#include "HMM.h"
#include "DNAUtil.h"



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
	GHMM_Model* mdl = alloc<GHMM_Model>();
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


mvec<const char*> HMM::GetGenicCodons()
{
	static const char* genicCodons[] = {
		"AAA", "AAC", "AAT", "AAG", "ACA", "ACC", "ACT", "ACG", "ATA", "ATC", "ATT", "ATG",
		"AGA", "AGC", "AGT", "AGG", "CAA", "CAC", "CAT", "CAG", "CCA", "CCC", "CCT", "CCG",
		"CTA", "CTC", "CTT", "CTG", "CGA", "CGC", "CGT", "CGG", "TAC", "TAT", "TCA", "TCC",
		"TCT", "TCG", "TTA", "TTC", "TTG", "TTT", "TGC", "TGT", "TGG", "GAA", "GAC", "GAT",
		"GAG", "GCA", "GCC", "GCT", "GCG", "GTA", "GTC", "GTT", "GTG", "GGA", "GGC", "GGT", "GGG"
	};
	return mvec<const char*>(genicCodons, sizeof(genicCodons)/sizeof(char*));
}

mvec<const char*> HMM::GetStartCodons()
{
	static const char* startCodons[] = {"ATG", "GTG", "TTG"};
	return mvec<const char*>(startCodons, sizeof(startCodons)/sizeof(char*));
}

mvec<const char*> HMM::GetStopCodons()
{
	static const char *stopCodons[] = {"TAA", "TGA", "TAG"};
	return mvec<const char*>(stopCodons, sizeof(stopCodons)/sizeof(char*));
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

void HMM::AddStopCodons()
{
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

}

void HMM::AddStartCodons()
{

}

void HMM::TestModel()
{

	fprintf(stdout,"transition matrix:\n");
	model_A_print(stdout,ghmm_mdl,""," ","\n");
	fprintf(stdout,"observation symbol matrix:\n");
	model_B_print(stdout,ghmm_mdl,""," ","\n");

	sequence_t *seq = model_generate_sequences(ghmm_mdl,0,20,1,100);
	sequence_print(stdout, seq);
	sequence_free(&seq);
}
// 
//         % Returns start codon statistics for given sequence
//         % Input:
//         %       seq                    - sequence with annotated genes
//         %       traditionalNucleotides - boolean switch to return only
//         %                                genes with traditional nucleotides
//         %       traditionalCodons      - boolean switch to return only
//         %                                genes with traditional start and
//         %                                stop codons
//         %       length3                - boolean switch to return only
//         %                                genes with length devisible by 3
//         % Output:
//         %       startC                 - array with codon statistics for
//         %                                start codons
//         %       stopC                  - array with codon statistics for
//         %                                stop codons
//         function [startC, stopC] = get_start_stop_statistics(seq, traditionalNucleotides, length3, traditionalCodons)
//             startC = zeros([16 16 16]);
//             stopC = zeros([16 16 16]);
//             n = length(seq.gene);
//             for i = 1:n
//                 st = seq.gene(i).Indices(1);
//                 fn = seq.gene(i).Indices(2);
//                 len = abs(fn - st) + 1;
//                 if (~length3 || mod(len, 3) == 0)
//                     rev = false;
//                     if (st < fn)
//                         start = seq.Sequence(st:st + 2);
//                         stop = seq.Sequence(fn - 2:fn);
//                         gene = seq.Sequence(st:fn);
//                     else
//                         start = seq.Sequence(st - 2:st);
//                         stop = seq.Sequence(fn:fn + 2);
//                         gene = seqrcomplement(seq.Sequence(fn:st));
//                         rev = true;
//                     end
//                     
//                     if (~traditionalNucleotides || sum(nt2int(gene) > 4) == 0)
//                         if (rev == true)
//                             start = seqrcomplement(start);
//                             stop = seqrcomplement(stop);
//                         end
//                         if (~traditionalCodons || (HMM.is_start(start) && HMM.is_stop(stop)))
//                             startC(nt2int(start(1)), nt2int(start(2)), nt2int(start(3))) = startC(nt2int(start(1)), nt2int(start(2)), nt2int(start(3))) + 1;
//                             stopC(nt2int(stop(1)), nt2int(stop(2)), nt2int(stop(3))) = stopC(nt2int(stop(1)), nt2int(stop(2)), nt2int(stop(3))) + 1;
//                         end
//                     end
//                 end
//             end
//         end

// 
        //         % Adds stop codons to the model
        //         function [] = add_stop(obj)
        //             id1(1) = obj.add_state('stop_TAA_TGA_1', HMM.get_emission_prob_for_nucleotide('T'));
        //             count_TAA = obj.stop_codon_stats(nt2int('T'), nt2int('A'), nt2int('A'));
        //             count_TGA = obj.stop_codon_stats(nt2int('T'), nt2int('G'), nt2int('A'));
        //             count_both = count_TAA + count_TGA;
        //             id1(2) = obj.add_state('stop_TAA_TGA_2', (HMM.get_emission_prob_for_nucleotide('A') * count_TAA + HMM.get_emission_prob_for_nucleotide('G') * count_TGA) ./ count_both);
        //             id1(3) = obj.add_state('stop_TAA_TGA_3', HMM.get_emission_prob_for_nucleotide('A'));
        //             obj.add_edge(id1(1), id1(2), 1);
        //             obj.add_edge(id1(2), id1(3), 1);
        //             
        //             id2(1) = obj.add_state('stop_TAG_1', HMM.get_emission_prob_for_nucleotide('T'));
        //             id2(2) = obj.add_state('stop_TAG_2', HMM.get_emission_prob_for_nucleotide('A'));
        //             id2(3) = obj.add_state('stop_TAG_3', HMM.get_emission_prob_for_nucleotide('G'));
        //             
        //             obj.add_edge(id2(1), id2(2), 1);
        //             obj.add_edge(id2(2), id2(3), 1);
        // 
        // %              % Single entry point to stop codons!
        // %              id_main = obj.add_state('stop_T', HMM.get_emission_prob_for_nucleotide('T'));
        // %               
        // %              count_TAA = obj.stop_codon_stats(nt2int('T'), nt2int('A'), nt2int('A'));
        // %              count_TGA = obj.stop_codon_stats(nt2int('T'), nt2int('G'), nt2int('A'));
        // %              count_both = count_TAA + count_TGA;
        // %              id1(1) = obj.add_state('stop_TAA_TGA_2', (HMM.get_emission_prob_for_nucleotide('A') * count_TAA + HMM.get_emission_prob_for_nucleotide('G') * count_TGA) ./ count_both);
        // %              id1(2) = obj.add_state('stop_TAA_TGA_3', HMM.get_emission_prob_for_nucleotide('A'));
        // %              obj.add_edge(id1(1), id1(2), 1);
        // %              
        // %              count_TAG = obj.stop_codon_stats(nt2int('T'), nt2int('A'), nt2int('G'));
        // %              id2(1) = obj.add_state('stop_TAG_2', HMM.get_emission_prob_for_nucleotide('A'));
        // %              id2(2) = obj.add_state('stop_TAG_3', HMM.get_emission_prob_for_nucleotide('G'));
        // %              obj.add_edge(id2(1), id2(2), 1);
        // %               
        // %              % Now to link single entry point to 2 stop codon models
        // %              count_all = count_both + count_TAG;
        // %              obj.add_edge(id_main, id1(1), count_both / count_all);
        // %              obj.add_edge(id_main, id2(1), count_TAG / count_all);
        //         end
        //         
        //         % Adds start codons to the model
        //         function [] = add_start(obj)
        //             count_ATG = obj.start_codon_stats(nt2int('A'), nt2int('T'), nt2int('G'));
        //             count_GTG = obj.start_codon_stats(nt2int('G'), nt2int('T'), nt2int('G'));
        //             count_TTG = obj.start_codon_stats(nt2int('T'), nt2int('T'), nt2int('G'));
        //             count_all = count_ATG + count_GTG + count_TTG;
        //             id(1) = obj.add_state('start_codons_AGT', (HMM.get_emission_prob_for_nucleotide('A') * count_ATG + HMM.get_emission_prob_for_nucleotide('T') * count_TTG + HMM.get_emission_prob_for_nucleotide('G') * count_GTG) ./ count_all);
        //             id(2) = obj.add_state('start_codons_T', HMM.get_emission_prob_for_nucleotide('T'));
        //             id(3) = obj.add_state('start_codons_G', HMM.get_emission_prob_for_nucleotide('G'));
        //             obj.add_edge(id(1), id(2), 1);
        //             obj.add_edge(id(2), id(3), 1);
//         end
