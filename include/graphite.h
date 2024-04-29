#ifndef __IMPUTEVCF
#define __IMPUTEVCF

#include <map>
#include <set>
#include "VCFHeteroHomo.h"
#include "VCFFillable.h"
#include "Pedigree.h"
#include "Map.h"

class Option;
class SampleManager;
class VCFFamily;
class VCFHeteroHomo;


//////////////////// Materials ////////////////////

class Materials {
	const std::string	path_map;
	const Map	*geno_map;
	const std::vector<const Map *>	chr_maps;
	
public:
	Materials(const std::string& path, const Map *m);
	~Materials();
	
	const std::string& get_path_map() const { return path_map; }
	const Map& get_map() const { return *geno_map; }
	const Map *get_chr_map(int i) const;
	double total_cM() const;
	void display_map_info() const;
	
public:
	static Materials *create(const Option *option);
};


//////////////////// process ////////////////////

void display_chromosome_info(const VCFSmall *orig_vcf);
VCFSmall *impute_vcf_chr(const VCFSmall *vcf, SampleManager *sample_man,
								const Map& geno_map, const Option *option);
void print_info(const Option *option);
void impute_all(VCFHuge *vcf, const Materials *materials, const Option *option);
void impute_progenies(VCFHuge *vcf, const Materials *materials,
												const Option *option);
void impute_VCF(const Option *option);
#endif
