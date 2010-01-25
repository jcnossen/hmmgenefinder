#pragma once

class CfgList;

using std::vector;
using std::string;

class Feature
{
public:
	Feature(CfgList* d);
	virtual Feature* Clone();

	bool Complementary() { return indices[0] > indices[1]; }

	mvec<int> indices;
	std::string gene, note, func, locusTag, product;

	enum Type {
		Type_tRNA, Type_rRNA, Type_MiscRNA, Type_MiscFeature, Type_Protein
	} type;
};

// Make mvec.clone() use the proper cloning method
template<> class cloning_device<Feature> {
public:
	static Feature* clone(Feature *src) { return src->Clone(); }
};

class FeatureProtein : public Feature {
public:
	FeatureProtein(CfgList* d);
	Feature* Clone();

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
	void PrintGenes();

	// returns dna or reverse complement dna
	std::string GetGeneDNA(Feature* f);

	mvec< mvec<int>* > GetGenicDNA();
	mvec< mvec<int>* > GetIntergenicDNA();

	// all features combined
	mvec<Feature*> genes;
	
	string sequence;
	string filename;

	string name, topology, molecularType;
	string modificationDate, source, definition;

	mvec<Genome*> Split(float wantedRatio=0.5f, int impTh=5);

	Genome* GetSubset(int first, int last);
	Genome* GetSubsetByGeneIndex(int first, int last);

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
