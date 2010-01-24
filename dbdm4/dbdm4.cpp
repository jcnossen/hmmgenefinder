// dbdm4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "Genome.h"

#include "math.h"
#include "HMM.h"
#include "DNAUtil.h"



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

	mvec<int> stateSeq;
	mvec<int> testSeq = hmm->GenerateSequence(1000, &stateSeq);
	d_trace("Test sequence (length %d)\n", testSeq.size());
//  	for (int i=0;i<testSeq.size();i++)
//  		d_trace("\t[%d]=%d\n", i,testSeq[i]);

	hmm->BaumWelch(testSeq);

	mvec<int> viterbiPath = hmm->ViterbiPath(testSeq);

 	d_trace("Viterbi path (%d): \n", viterbiPath.size());
 	for(int i=0;i<viterbiPath.size();i++)
 		d_trace("State: %s\n", hmm->states[viterbiPath[i]]->name.c_str());

	delete hmm;

	getc(stdin);
}




int main(int argc, char* argv[])
{
	try {
		Genome genome ("../data/AE005174.gd");
		genome.sequence = dna::RandomizeUnknownNT(genome.sequence);

		// select a small piece of genome
		Genome* piece = genome.GetSubsetByGeneIndex(0, 150);

		SimpleGenomicTest test;

		mvec<Genome*> tt = piece->Split();
		Genome *train = tt(1), *test_genome = tt(2);
		d_trace("Train genome: \n");
		train->PrintInfo();
		d_trace("Test genome: \n");
		test_genome->PrintInfo();

		test.Build();
		test.Train(train);


		delete genicHMM;
		DeleteAll(tt);
	}
	catch (ContentException e) {
		d_trace("Exception: %s\n", e.what());
	}
	return 0;
}

