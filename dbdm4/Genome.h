#pragma once

class CfgList;

using std::vector;
using std::string;

class Feature
{
public:
	Feature(CfgList* d);

	mvec<int> indices;
	std::string gene, note, func, locusTag, product;
};

class FeatureProtein : public Feature {
public:
	FeatureProtein(CfgList* d);

	string proteinID, translation; 
	int codonStart;
};

class Genome
{
public:
	Genome() {}
	Genome(std::string file);
	~Genome();

	void PrintInfo();

	mvec<FeatureProtein*> proteins;
	mvec<Feature*> tRNA, rRNA, misc_RNA, misc_feature;

	// all features combined
	mvec<Feature*> genes;
	
	string sequence;
	string filename;

	string name, topology, molecularType;
	string genBankDivision, modificationDate, source, definition;

	mvec<Genome*> Split(float wantedRatio=0.5f, int impTh=5);

	Genome* GetSubset(int first, int last);

private:

	// Returns a boolean value - tells us if we can cut the genome at given
	// position.
	bool CanCut(int position);


	// Returns all genes that start and end in given range
	// src   - genes to select from
	// s,f - start and end of allowed nucleotide range
	template<typename T> static mvec<T*> GetGenesInNucleotideRange(mvec<T*> src, int s, int f);

	template<typename T> mvec<T*> LoadFeatureList(CfgList* l);


// 	% Returns number of genes that start and end in given range
	int CountGenes(int start, int finish);

};
