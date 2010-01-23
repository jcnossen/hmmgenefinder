// dbdm4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "Genome.h"

#include "math.h"
#include "HMM.h"


int random() { return rand(); }
void srandom(unsigned int s) { srand(s); }


void TestHMM()
{
	HMM *hmm = new HMM();

	float f[] = { 0.2f, 0.8f, 0.0f, 1.0f };
	HMMState* s = hmm->AddState("coinflip", mvec<float>(f, 4));

	s->AddEdge(s, 0.8f);

	hmm->NormalizeProbabilities();
	hmm->BuildModel();

	hmm->TestModel();

	getc(stdin);
}

int main(int argc, char* argv[])
{
	try {
		TestHMM();
/*		Genome genome ("../data/AE005174.gd");

		genome.PrintInfo();

		mvec<Genome*> tt = genome.Split();
		Genome *train = tt(1), *test = tt(2);

		d_trace("Train genome: \n");
		train->PrintInfo();
		d_trace("Test genome: \n");
		test->PrintInfo();

		DeleteAll(tt);*/
	}
	catch (ContentException e) {
		d_trace("Exception: %s\n", e.what());
	}
	return 0;
}

