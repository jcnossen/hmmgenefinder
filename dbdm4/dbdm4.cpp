// dbdm4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "Genome.h"

#include "math.h"
#include "HMM.h"
#include "DNAUtil.h"


std::string GetGenicDNA(Genome* genome) 
{
	std::string seq;
	// gather codon stats for genic regions
	for (int i=0;i<genome->genes.size();i++) {
		Feature* f = genome->genes[i];

		seq += genome->GetGeneDNA(f);
	}
	return seq;
}

HMM* CreateGenicHMM( Genome* genome )
{
	mvec<int> cc(64, 0);

	// gather codon stats for genic regions
	for (int i=0;i<genome->genes.size();i++) {
		Feature* f = genome->genes[i];

		//	d_trace("CC: [%d,%d] len: %d: %s\n", f->indices[0], f->indices[1], abs(f->indices[1]-f->indices[0]), f->gene.c_str());

		std::string seq = genome->GetGeneDNA(f);
		mvec<int> seqCC = dna::CodonCount(seq);
		cc += seqCC;
	}

	//	dna::ListCodonCounts(cc);

	HMMState* center = 0;

	HMM* hmm =  new HMM();
	center = hmm->AddState("center");

	for (int i=0;i<64;i++) {
		dna::Codon c(i);
		HMMState* prev = 0;
		mvec<float> emit (4);
		for (int j=0;j<3;j++) {
			int nt = dna::nt2int(c.nt[j]);

			for (int x=0;x<4;x++)
				emit[x] = x==nt?1.0f:0.0f;

			HMMState* n = hmm->AddState(std::string(c.nt) + "_" + c.nt[j], emit);
			if (prev)
				prev->AddEdge(n, 1.0f);
			else
				center->AddEdge(n, cc[i]);
			prev = n;
		}
		//connect last nucleotide back to center
		prev->AddEdge(center, 1.0f);
	}

	return hmm;
}



void TestHMM()
{
	HMM *hmm = new HMM();

	float s_e[] = { 1,1,1,1 }; // random emissions
	float c1_e[] = { 0,1,0,0 }; // always 1
	float c3_e[] = { 0,0,0,1 }; // always 3

	HMMState* s = hmm->AddState("rnd", mvec<float>(s_e, 4));
	HMMState* c1 = hmm->AddState("c1", mvec<float>(c1_e, 4));
	HMMState* c3 = hmm->AddState("c3", mvec<float>(c3_e, 4));

	// Edges to self
	s->AddEdge(s, 0.9f); // 0.9 chance it stays random
	c1->AddEdge(c1, 0.9f);
	c3->AddEdge(c3, 0.9f);

	// Edges to others
	c1->AddEdge(c3, 0.1f);
	c3->AddEdge(s, 0.1f);
	s->AddEdge(c1, 0.1f);

	hmm->BuildModel();

	//hmm->TestModel();

	mvec<int> testSeq = hmm->GenerateSequence(1000);
	d_trace("Test sequence (length %d)\n", testSeq.size());
//  	for (int i=0;i<testSeq.size();i++)
//  		d_trace("\t[%d]=%d\n", i,testSeq[i]);

	hmm->BaumWelch(testSeq);

	mvec<int> viterbiPath = hmm->ViterbiPath(testSeq);

// 	d_trace("Viterbi path (%d): \n", viterbiPath.size());
// 	for(int i=0;i<viterbiPath.size();i++)
// 		d_trace("State: %s\n", hmm->states[viterbiPath[i]]->name.c_str());

	delete hmm;

	getc(stdin);
}




int main(int argc, char* argv[])
{
	try {
		Genome genome ("../data/AE005174.gd");
		genome.sequence = dna::RandomizeUnknownNT(genome.sequence);

		// select a small piece of genome
		Genome* piece = genome.GetSubsetByGeneIndex(0, 200);

		HMM* genicHMM = CreateGenicHMM(piece);

		piece->PrintInfo();

		mvec<Genome*> tt = piece->Split();
		Genome *train = tt(1), *test = tt(2);

		genicHMM->BuildModel();

		std::string train_genic = GetGenicDNA(train);
		genicHMM->BaumWelch(dna::nt2int(train_genic));

		d_trace("Train genome: \n");
		train->PrintInfo();
		d_trace("Test genome: \n");
		test->PrintInfo();

		DeleteAll(tt);
	}
	catch (ContentException e) {
		d_trace("Exception: %s\n", e.what());
	}
	return 0;
}

