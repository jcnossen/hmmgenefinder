
#include "stdafx.h"

#include "DNAUtil.h"

namespace dna {

// change all unknown nucleotides to a random nucleotide
std::string RandomizeUnknownNT(std::string dna) 
{
	for (int i=0;i<dna.size();i++) {
		char c = dna[i];
		if (c != 'a' && c != 'c' && c != 'g' && c!='t')
			dna[i] = ("acgt") [rand() & 3];
	}
	return dna;
}

int nt2int( char v )
{
	 if (v == 'a') return 0;
	 if (v == 'c') return 1;
	 if (v == 'g') return 2;
	 if (v == 't' || v=='u') return 3;
	 throw ContentException(SPrintf("Invalid char passed to nt2int:%c",v));
}

mvec<int> nt2int( const std::string& v )
{
	mvec<int> r(v.size());
	for (int x=0;x<v.size();x++)
		r[x]=nt2int(v[x]);
	return r;
}



std::string SeqComplement( std::string dna )
{
	char compl_map[5] = {'t', 'g', 'c', 'a' }; // mapped to return values of nt2int
	std::string compl;

	for (int i=0;i<dna.size();i++)
		compl += compl_map[nt2int(dna[i])];
	return compl;
}

std::string SeqReverseComplement( std::string dna )
{
	return SeqComplement(std::string (dna.rbegin(), dna.rend()));
}

char int2nt( int i )
{
	const char *c="acgt";
	return c[i];
}

std::string int2nt( const mvec<int>& v )
{
	std::string r;
	r.resize(v.size());
	for(int i=0;i<v.size();i++)
		r[i]=int2nt(v[i]);
	return r;
}


Codon::Codon( const std::string& v )
{
	assert(v.size()==3);
	nt[0]=v[0];
	nt[1]=v[1];
	nt[2]=v[2];
	nt[3]=0;
}

Codon::Codon( const char* v )
{
	for (int i=0;i<3;i++)
		nt[i]=v[i];
	nt[3] = 0;
}


static const char* startc[] = {"atg", "gtg", "ttg"};
static const char *stopc[] = {"taa", "tga", "tag"};


CodonType GetCodonType( int id )
{
	static int codonMap[64] = {-1};

	if (codonMap[0]<0) {
		for(int i=0;i<64;i++)
			codonMap[i] = Codon_Genic;

		for (int i=0;i<sizeof(startc)/sizeof(char*);i++) {
			codonMap[ Codon(startc[i]).GetID() ] = Codon_Start;
		}
		for (int i=0;i<sizeof(stopc)/sizeof(char*);i++) {
			codonMap[ Codon(stopc[i]).GetID() ] = Codon_Stop;
		}
		// the rest of it is genic
	}

	return (CodonType)codonMap[id];
}

void ListCodonCounts(mvec<int> cc) {
	for (int i=0;i<cc.size();i++){
		d_trace("%s: %d\n", Codon(i).nt, cc[i]);
	}
}

mvec<int> CodonCount( std::string seq )
{
	mvec<int> r;
	r.resize(64);
	
	int p = 0;
	while (p < seq.length() - 2) {
		Codon c(&seq[p]);
//		d_trace("Codon:[%d] %s\n", c.GetID(),c.nt); 
		r[c.GetID()] ++;
		p += 3;
	}
	return r;
}

};
