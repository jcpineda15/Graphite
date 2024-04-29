#ifndef __SMALLFAMILY
#define __SMALLFAMILY

#include <vector>
#include "../include/VCF.h"

class VCFSmall;
class VCFFamily;
class VCFFillable;
class SampleManager;
class KnownFamily;
class Map;
class Option;

namespace SmallFamily {
	VCFRecord *merge_progeny_records(std::vector<VCFFillable *>& vcfs,
									std::size_t i, const STRVEC& samples);
	
	VCFSmall *impute_vcf_by_parents_core(
						const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
						const std::vector<const KnownFamily *>& families,
						const Map& geno_map, const Option *option);
	VCFSmall *impute_vcf_by_parents(
						const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
						SampleManager *sample_man,
						const Map& geno_map, const Option *option);
	
	VCFSmall *impute_vcf_by_parent_core(
						const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
						const std::vector<const KnownFamily *>& families,
						const Map& geno_map,
						SampleManager *sample_man, const Option *option);
	VCFSmall *impute_vcf_by_parent(const VCFSmall *orig_vcf,
							const VCFSmall *merged_vcf, const Map& geno_map,
							SampleManager *sample_man, const Option *option);
	
	VCFSmall *impute_one_parent_vcf_core(const VCFSmall *orig_vcf,
							const VCFSmall *merged_vcf,
							const std::vector<const KnownFamily *>& families,
							const Map& geno_map,
							SampleManager *sample_man, int num_threads);
	VCFSmall *impute_one_parent_vcf(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads);
	
	VCFSmall *impute_small_family_VCFs(const VCFSmall *orig_vcf,
										VCFSmall *merged_vcf,
										const Map& geno_map,
										SampleManager *sample_man,
										const Option *option);
	
	VCFSmall *impute_vcf_by_progenies_core(const VCFSmall *orig_vcf,
							const VCFSmall *merged_vcf,
							const std::vector<const KnownFamily *>& families,
							const Map& geno_map,
							SampleManager *sample_man, int num_threads);
	VCFSmall *impute_vcf_by_progenies(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads);
	
	VCFSmall *impute_iolated_samples(
					const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
					SampleManager *sample_man, const STRVEC& samples,
					const Map& gmap, bool modify_genotypes, int num_threads);
}
#endif
