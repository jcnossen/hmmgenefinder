
#include "stdafx.h"

#include "DNAUtil.h"



int nt2int( char v )
{
// A => 1 , C => 2, G => 3,  T(U) => 4.
	 if (v == 'A') return 1;
	 if (v == 'C') return 2;
	 if (v == 'G') return 3;
	 if (v == 'T' || v=='U') return 4;
	 return 0;
}

mvec<int> nt2int( std::string v )
{
	mvec<int> r(v.size());
	for (int x=0;x<v.size();x++)
		r[x]=nt2int(v[x]);
	return r;
}

std::string seqcomplement( std::string dna )
{
	char compl_map[5] = {' ', 'T', 'G', 'C', 'A' }; // mapped to return values of nt2int
	std::string compl;

	for (int i=0;i<dna.size();i++)
		compl += compl_map[nt2int(dna[i])];
	return compl;
}

std::string segrcomplement( std::string dna )
{
	return seqcomplement(std::string (dna.rend(), dna.rbegin()));
}

