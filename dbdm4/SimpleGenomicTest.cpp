#include "stdafx.h"
#include "DNAUtil.h"
#include "SimpleGenomicTest.h"

#include "HMM.h"
#include "Genome.h"

HMM* CreateIntergenicHMM( Genome* genome )
{
	// create an intergenic sequence
	std::string intergenic;

	for (int i=0;i<genome->genes.size();i++) {
		Feature* f = genome->genes[i];
	}

	return 0;
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

