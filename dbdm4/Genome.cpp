
#include "stdafx.h"

#include "Genome.h"
#include "CfgParser.h"

int round(float f) { return (int)(f+0.5f); }

static std::string fz(const char* s) { return s?s:""; }

Genome::Genome(std::string file)
{
	filename = file;
	CfgList* d = CfgList::LoadFile(file);

	// read sequence data
	const char* seqfile = d->GetLiteral("SequenceFile");
	std::string seqPath = GetDirPath(file) + seqfile;
	sequence = ReadTextFile(seqPath);

	// load global properties

	name = d->GetLiteral("LocusName");
	topology = d->GetLiteral("LocusTopology");
	genBankDivision = d->GetLiteral("LocusGenBankDivision");
	modificationDate = d->GetLiteral("LocusModificationDate");
	definition = d->GetLiteral("Definition");
	source = d->GetLiteral("Source");

	// load features
	CfgList* features = d->GetList("features");
	mvec<FeatureProtein*> proteins = LoadFeatureList<FeatureProtein>(features->GetList("CDS"));
	members_ptr(proteins, &FeatureProtein::type) |= FeatureProtein::Type_Protein;
	
	mvec<Feature*> tRNA = LoadFeatureList<Feature>(features->GetList("tRNA"));
	members_ptr(tRNA, &Feature::type).ptr_assign(Feature::Type_tRNA);

	mvec<Feature*> rRNA = LoadFeatureList<Feature>(features->GetList("rRNA"));
	members_ptr(rRNA, &Feature::type).ptr_assign(Feature::Type_rRNA);

	mvec<Feature*> misc_RNA = LoadFeatureList<Feature>(features->GetList("misc_RNA"));
	members_ptr(misc_RNA, &Feature::type).ptr_assign(Feature::Type_MiscRNA);

	mvec<Feature*> misc_feature = LoadFeatureList<Feature>(features->GetList("misc_feature"));
	members_ptr(misc_feature, &Feature::type).ptr_assign(Feature::Type_MiscFeature);

	genes = tRNA & rRNA & misc_RNA & misc_feature & proteins;

	delete d;
}


template<typename T> 
mvec<T*> Genome::LoadFeatureList( CfgList* l )
{
	mvec<T*> features;
	features.reserve(l->childs.size());

	for(CfgList::iterator li = l->begin(); li != l->end(); ++li) {
		features.push_back(new T( (CfgList*)(*li)->value ));
	}

	return features;
}

Genome::~Genome()
{
	DeleteAll(genes); // genes contains the objects from all other containers as well
}

void Genome::PrintInfo()
{
	d_trace("Genome %s contains: \n", name.c_str());
	d_trace("\t%d nucleotides\n", sequence.size());

	d_trace("\t%d proteins\n", sum(members(genes, &Feature::type) == Feature::Type_Protein));
	d_trace("\t%d tRNA parts\n", sum(members(genes, &Feature::type) == Feature::Type_tRNA));
	d_trace("\t%d rRNA parts\n", sum(members(genes, &Feature::type) == Feature::Type_rRNA));
	d_trace("\t%d misc. RNA\n", sum(members(genes, &Feature::type) == Feature::Type_MiscRNA));
	d_trace("\t%d misc. features\n", sum(members(genes, &Feature::type) == Feature::Type_MiscFeature));
}


// Returns a boolean value - tells us if we can cut the genome at given
// position.
bool Genome::CanCut(int location) {
	for(int i=1;i<=genes.size();i++) {
		Feature* g = genes(i);
		if (min(g->indices) <= location && location <= max(g->indices))
			return false;
	}
	return true;
}

// 	% Returns all genes that start and end in given range
// 		% g   - genes to select from
// 		% s,f - start and end of allowed nucleotide range
template<typename T>
mvec<T*> Genome::GetGenesInNucleotideRange(mvec<T*> src, int s, int f) 
{
	mvec<T*> genes;
	for(int gagI = 1; gagI <= src.size(); ++gagI) {
		if (min(src(gagI)->indices) >= s && max(src(gagI)->indices) <= f)
			genes.push_back(src(gagI));
	}
	return genes;
}


mvec<Genome*> Genome::Split( float wanted_ratio, int impTh )
{
	/*
	% For each gene go and try to divide the genome after each gene, make
	% sure we do not cut any genes in the middle and select the ratio
	% closest to 0.5 */
	float best_ratio = FLT_MAX;
	int best_position = 0;
	int n = genes.size();

	int last_impI = 0;
	int last_impJ = 0;
	int i = (int)(n * wanted_ratio + 0.5f) - 1;
	int j = (int)(n * wanted_ratio + 0.5f);
	while ((i > 0 && last_impI < impTh) || (j < n && last_impJ < impTh))
	{
		if (i > 0 && last_impI < impTh) {
			Feature* cur = genes(i);
			Feature* next = genes(i + 1);
			int cur_end = max(cur->indices);
			int next_start = min(next->indices);
			int middle = round(0.5 * (cur_end+next_start));
			if (CanCut(middle)) {
				float ratio = CountGenes(1, middle) / (float)n;
				if (fabsf(ratio - wanted_ratio) < fabsf(best_ratio - wanted_ratio)) {
					best_ratio = ratio;
					d_trace("[+] (%d) New best ratio attained - %f\n", i, best_ratio);
					best_position = middle;
					last_impI = 0;
				} else {
					last_impI = last_impI + 1;
				}
			}
		}
		i--;

		if (j < n && last_impJ < impTh) {
			Feature* cur = genes(j);
			Feature* next = genes(j + 1);
			int cur_end = max(cur->indices);
			int next_start = min(next->indices);
			int middle = round(0.5 * (cur_end+next_start));
			if (CanCut(middle)) {
				float ratio = CountGenes(1, middle) / (float) n;
				if (fabsf(ratio - wanted_ratio) < fabsf(best_ratio - wanted_ratio)) {
					best_ratio = ratio;
					d_trace("[+] (%d) New best ratio attained - %f\n", i, best_ratio);
					best_position = middle;
					last_impJ = 0;
				} else {
					last_impJ = last_impJ + 1;
				}
			}
		}
		j++;
	}
// 		% BTW, this works only coz the genes are sorted in incresing order of
// 		% their lower index (lower != first)

	d_trace("[i] Cutting sequence at %d\n", best_position); 
	mvec<Genome*> r;
	r.push_back(GetSubset(1, best_position)); // train
//	train.Sequence = g.Sequence(1:best_position);
	//train.gene = get_all_genes(f, 1, best_position);
	r.push_back(GetSubset(best_position + 1, sequence.size()));
// 	test.Sequence = g.Sequence(best_position + 1:seq_length);
// 	test.gene = get_all_genes(f, best_position + 1, seq_length);
// 	test.gene = shift_genes(test.gene, best_position);
	return r;
}

Genome* Genome::GetSubset( int s, int f)
{
	Genome* g = new Genome();
	g->sequence = sequence.substr(s, f-s);
	g->definition = definition;
	g->source = source;
	g->genBankDivision = genBankDivision;
	g->topology = topology;
	g->name = name;
	g->molecularType = molecularType;
	g->filename = filename;
	g->modificationDate = modificationDate;

	g->genes = GetGenesInNucleotideRange<Feature>(genes, s,f).clone();

	// shift genes
	for(int i=1;i<=g->genes.size();++i)
		g->genes(i)->indices = g->genes(i)->indices - (s-1); // FIXED? it should be the nucleotides before best_position, which is best_position-1

	return g;
}

int Genome::CountGenes( int start, int finish )
{
	int cnt = 0;
	for(int cgI=1;cgI<=genes.size();++cgI)
		if (min(genes(cgI)->indices) >= start && max(genes(cgI)->indices) <= finish)
			cnt = cnt + 1;
	return cnt;
}

FeatureProtein::FeatureProtein( CfgList* d ) : Feature(d)
{
	codonStart = d->GetInt("codon_start");
	proteinID = d->GetLiteral("protein_id");
	translation = d->GetLiteral("translation");

	type = Feature::Type_Protein;
}

Feature* FeatureProtein::Clone()
{
	return new FeatureProtein(*this);
}

Feature::Feature( CfgList* d )
{
	gene = fz(d->GetLiteral("gene"));
	func = fz(d->GetLiteral("function"));
	locusTag = fz(d->GetLiteral("locus_tag"));
	note = fz(d->GetLiteral("note"));
	product = fz(d->GetLiteral("product"));

	CfgList* indexList = d->GetList("Indices");
	indices.push_back((int)((CfgNumeric*) indexList->childs[0]->value)->value);
	indices.push_back((int)((CfgNumeric*) indexList->childs[1]->value)->value);
}

Feature* Feature::Clone()
{
	return new Feature(*this);
}
