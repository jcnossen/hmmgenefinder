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

	double s_e[] = { 1,1,1,1 }; // random emissions
	double c1_e[] = { 0,1,0,0 }; // always 1
	double c3_e[] = { 0,0,0,1 }; // always 3

	HMMState* s = hmm->AddState("rnd", mvec<double>(s_e, 4));
	HMMState* c1 = hmm->AddState("c1", mvec<double>(c1_e, 4));
	HMMState* c3 = hmm->AddState("c3", mvec<double>(c3_e, 4));

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

	mvec<int> viterbiPath = hmm->ViterbiPath(testSeq);

 	d_trace("Viterbi path (%d): \n", viterbiPath.size());
 	for(int i=0;i<viterbiPath.size();i++)
 		d_trace("State: %s\n", hmm->states[viterbiPath[i]]->name.c_str());

	delete hmm;

	getc(stdin);
}


mvec< mvec<int>* > load_sequence_set(string file)
{
	std::string text = ReadTextFile(file);
	mvec< mvec<int>* > results; 

	int pos = 0;
	while (pos < text.length()) {
		int start = pos;
		while (isalpha(text[pos]) && pos < text.length()) 
			pos ++;

		std::string dnastr = text.substr(start, pos-start-1);
		dnastr = dna::RandomizeUnknownNT(dnastr);
		results.push_back( new mvec<int>(dna::nt2int(dnastr)) );
		
		while (!isalpha(text[pos]) && pos < text.length())
			pos ++;
	}
	return results;
}

HMM* CreateGenicHMM( Genome* genome );

void TestGenomic()
{
	Genome genome ("../bin/AE005174.gd");
	genome.sequence = dna::RandomizeUnknownNT(genome.sequence);

	// select a small piece of genome
	Genome* piece = genome.GetSubsetByGeneIndex(0, 50);
	mvec<Genome*> tt = piece->Split();
	Genome *train = tt(1), *test_genome = tt(2);
	d_trace("Train genome: \n");
	train->PrintInfo();
	d_trace("Test genome: \n");
	test_genome->PrintInfo();

	HMM *hmm = CreateGenicHMM(train);
	mvec< mvec<int>* > seqDNA = train->GetGenicDNA();
	hmm->BaumWelch(seqDNA);

	delete piece;
	DeleteAll(tt);
	DeleteAll(seqDNA);

}

void PrintHelp()
{
	d_trace("HMMUtil <options> <hmmfile> <sequencefile>\n"
		"Options:\n"
		"-viterbi\tReturn viterbi path of given sequences\n"
		"-bw\t\tTrain model using Baum Welch\n"
		);
}

int main(int argc, char* argv[])
{
	try {
		TestGenomic();

		

		return 0;
		int cmd = -1;

		if (argc < 4) {
			PrintHelp();
			return 0;
		}

		if (!STRCASECMP(argv[1], "-bw"))
			cmd = 1;
		else if (!STRCASECMP(argv[1], "-viterbi"))
			cmd = 2;
		else {
			PrintHelp();
			return 0;
		}

		std::string hmmFile = argv[2];
		std::string seqFile = argv[3];

		HMM* hmm = new HMM();
		hmm->ParseConfig(hmmFile);
		hmm->initial_state = hmm->FindState("start_codons_AGT");
		hmm->BuildModel();

		mvec< mvec<int> *> sequences = load_sequence_set(seqFile);

		if (cmd == 1) {
			hmm->BaumWelch(sequences);
		}
		if (cmd == 2) {
			for (int i=0;i<sequences.size();i++) {
				mvec<int> vp = hmm->ViterbiPath(*sequences[i]);
				printf("%d\n", vp.size());
				for (int j=0;j<vp.size();j++)
					printf("%d ", vp[j]);
				printf("\n");
			}
		}

		getc(stdin);

		DeleteAll(sequences);
	}
	catch (ContentException e) {
		d_trace("Exception: %s\n", e.what());
	}
	return 0;
}

