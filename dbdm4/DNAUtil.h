
#pragma once

namespace dna
{
	int nt2int(char v);
	char int2nt(int i);
	mvec<int> nt2int(const std::string& v);
	std::string int2nt(const mvec<int>& v);

	std::string SeqComplement(std::string dna);
	std::string SeqReverseComplement(std::string dna);
	std::string RandomizeUnknownNT(std::string dna);

	enum CodonType {
		Codon_Genic, Codon_Start, Codon_Stop
	};

	CodonType GetCodonType(int id);
	mvec<int> CodonCount(std::string seq);
	void ListCodonCounts(mvec<int> cc);

	struct Codon
	{
		Codon(const std::string& v);
		Codon(const char* v);

		Codon(int id) {
			nt[0] = int2nt( id & 3 );
			nt[1] = int2nt( (id / 4) & 3 );
			nt[2] = int2nt( (id / 16) & 3 );
			nt[3] = 0;
		}

		int GetID() {
			return nt2int(nt[0]) + nt2int(nt[1]) * 4 + nt2int(nt[2]) * 16;
		}

		CodonType GetType() { return GetCodonType(GetID()); }

		char nt[4]; // length 4, so it's usable as a C string
	};

};
