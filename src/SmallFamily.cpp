#include "../include/SmallFamily.h"
#include "../include/VCF.h"
#include "../include/VCFFillable.h"
#include "../include/VCFHeteroHomoPP.h"
#include "../include/VCFOneParentPhased.h"
#include "../include/VCFProgenyPhased.h"
#include "../include/VCFIsolated.h"
#include "../include/SampleManager.h"
#include "../include/common.h"

using namespace std;

VCFRecord *SmallFamily::merge_progeny_records(vector<VCFFillable *>& vcfs,
											size_t i, const STRVEC& samples) {
	const STRVEC&	v1 = vcfs.front()->get_records()[i]->get_v();
	STRVEC	v(v1.begin(), v1.begin() + 9);
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const STRVEC&	v2 = (*p)->get_records()[i]->get_v();
		v.insert(v.end(), v2.begin() + 11, v2.end());
	}
	return new VCFRecord(v, samples);
}

VCFSmall *SmallFamily::impute_vcf_by_parents_core(
						const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
						const vector<const KnownFamily *>& families,
						const Map& geno_map, const Option *option) {
	auto	vcfs = VCFHeteroHomoPP::impute_vcfs(orig_vcf, merged_vcf, families,
												geno_map, option->num_threads);
	STRVEC	samples;	// collect progenies
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		VCFFillable	*vcf = *p;
		STRVEC	ss = vcf->get_samples();
		samples.insert(samples.end(), ss.begin() + 2, ss.end());
	}
	
	// Give ownership of samples to new _vcf
	auto	new_header = orig_vcf->trim_header(samples);
	vector<VCFRecord *>	empty_records;
	auto	*new_vcf = new VCFSmall(new_header, samples, empty_records);
	const auto&	samples_ = new_vcf->get_samples();
	
	vector<VCFRecord *>	merged_records;
	for(size_t i = 0; i < vcfs.front()->size(); ++i) {
		new_vcf->add_record(merge_progeny_records(vcfs, i, samples_));
	}
	
	Common::delete_all(vcfs);
	
	return new_vcf;
}

VCFSmall *SmallFamily::impute_vcf_by_parents(
			const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
			SampleManager *sample_man,
			const Map& geno_map, const Option *option) {
	auto	families = sample_man->extract_small_families();
	if(families.empty())
		return NULL;
	
	auto	*new_imputed_vcf = impute_vcf_by_parents_core(orig_vcf,
												merged_vcf, families,
												geno_map, option);
	Common::delete_all(families);
	auto	*new_merged_vcf = VCFSmall::join(merged_vcf, new_imputed_vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(new_imputed_vcf->get_samples());
	delete new_imputed_vcf;
	return new_merged_vcf;
}

VCFSmall *SmallFamily::impute_vcf_by_parent_core(
					const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
					const vector<const KnownFamily *>& families,
					const Map& geno_map,
					SampleManager *sample_man, const Option *option) {
	// collect not phased parents
	STRVEC	samples;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		if(sample_man->is_imputed(family->get_mat()))
			samples.push_back(family->get_pat());
		else
			samples.push_back(family->get_mat());
	}
	
	// phase not phased parents
	VCFSmall	*parents_vcf = impute_iolated_samples(orig_vcf, merged_vcf,
														sample_man, samples,
														geno_map, true,
														option->num_threads);
	
	// merge vcfs
	VCFSmall	*new_merged_vcf = VCFSmall::join(merged_vcf, parents_vcf,
													orig_vcf->get_samples());
	
	// impute progenies
	VCFSmall	*new_vcf = impute_vcf_by_parents_core(orig_vcf, new_merged_vcf,
													families, geno_map, option);
	delete new_merged_vcf;
	
	// join
	VCFSmall	*vcf = VCFSmall::join(parents_vcf, new_vcf,
										orig_vcf->get_samples());
	delete parents_vcf;
	delete new_vcf;
	return vcf;
}

VCFSmall *SmallFamily::impute_vcf_by_parent(const VCFSmall *orig_vcf,
							const VCFSmall *merged_vcf, const Map& geno_map,
							SampleManager *sample_man, const Option *option) {
	auto	families = sample_man->extract_single_parent_phased_families();
	if(families.empty())
		return NULL;
	
	auto	*new_imputed_vcf = impute_vcf_by_parent_core(orig_vcf,
												merged_vcf, families, geno_map,
												sample_man, option);
	auto	*new_merged_vcf = VCFSmall::join(merged_vcf, new_imputed_vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(new_imputed_vcf->get_samples());
	delete new_imputed_vcf;
	Common::delete_all(families);
	return new_merged_vcf;
}

VCFSmall *SmallFamily::impute_one_parent_vcf_core(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const vector<const KnownFamily *>& families,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads) {
	STRVEC	references = sample_man->collect_large_family_parents();
	VCFSmall	*ref_vcf = merged_vcf->extract_samples(references);
	
	vector<VCFOneParentPhased *>	vcfs;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		const bool	is_mat_phased = sample_man->is_imputed(family->get_mat());
		auto	*vcf = VCFOneParentPhased::create(family->get_samples(),
												  is_mat_phased, merged_vcf,
												  orig_vcf, geno_map, ref_vcf);
		vcfs.push_back(vcf);
	}
	
	VCFSmall	*vcf = VCFOneParentPhased::impute_all(vcfs, num_threads);
	Common::delete_all(vcfs);
	delete ref_vcf;
	return vcf;
}

VCFSmall *SmallFamily::impute_one_parent_vcf(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads) {
	auto	families = sample_man->extract_phased_and_unknown_parents_family();
	if(families.empty())
		return NULL;
	
	auto	*new_imputed_vcf = impute_one_parent_vcf_core(orig_vcf,
										merged_vcf, families, geno_map,
										sample_man, num_threads);
	auto	*new_merged_vcf = VCFSmall::join(merged_vcf, new_imputed_vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(new_imputed_vcf->get_samples());
	delete new_imputed_vcf;
	Common::delete_all(families);
	return new_merged_vcf;
}

VCFSmall *SmallFamily::impute_vcf_by_progenies_core(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const vector<const KnownFamily *>& families,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads) {
	STRVEC	references = sample_man->collect_large_family_parents();
	VCFSmall	*ref_vcf = merged_vcf->extract_samples(references);
	
	vector<pair<const KnownFamily *, size_t>>	progeny_imputed_families;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		size_t	ppi;	// phased progeny index
		for(size_t i = 2; ; ++i) {
			if(sample_man->is_imputed(family->get_samples()[i])) {
				ppi = i;
				break;
			}
		}
		
		progeny_imputed_families.push_back(make_pair(family, ppi));
	}
	const auto	vcfs = VCFProgenyPhased::impute_all_by_progeny(
												orig_vcf, merged_vcf,
												progeny_imputed_families,
												geno_map, ref_vcf, num_threads);
	
	vector<const VCFSmallBase *>	vcfs2(vcfs.begin(), vcfs.end());
	VCFSmall	*new_vcf = VCFSmall::join(vcfs2, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	delete ref_vcf;
	return new_vcf;
}

// Impute families whose progenies have been imputed
VCFSmall *SmallFamily::impute_vcf_by_progenies(const VCFSmall *orig_vcf,
								  const VCFSmall *merged_vcf,
								  const Map& geno_map,
								  SampleManager *sample_man, int num_threads) {
	auto	families = sample_man->extract_progenies_phased_families();
	if(families.empty())
		return NULL;
	
	auto	*new_imputed_vcf = impute_vcf_by_progenies_core(orig_vcf,
													merged_vcf,
													families, geno_map,
													sample_man, num_threads);
	auto	*new_merged_vcf = VCFSmall::join(merged_vcf, new_imputed_vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(new_imputed_vcf->get_samples());
	delete new_imputed_vcf;
	Common::delete_all(families);
	return new_merged_vcf;
}

VCFSmall *SmallFamily::impute_small_family_VCFs(const VCFSmall *orig_vcf,
												VCFSmall *merged_vcf,
												const Map& geno_map,
												SampleManager *sample_man,
												const Option *option) {
	// Repeat until there are no more families to impute
	while(true) {
		// Impute families in which parents are imputed but children are few
		auto	new_merged_vcf1 = impute_vcf_by_parents(orig_vcf, merged_vcf,
														sample_man,
														geno_map, option);
		if(new_merged_vcf1 != NULL) {
			merged_vcf = new_merged_vcf1;
			continue;
		}
		
		// Impute families in which one parent is imputed
		auto	new_merged_vcf2 = impute_vcf_by_parent(orig_vcf, merged_vcf,
													   geno_map,
													   sample_man, option);
		if(new_merged_vcf2 != NULL) {
			merged_vcf = new_merged_vcf2;
			continue;
		}
		
		auto	*new_merged_vcf3 = impute_one_parent_vcf(orig_vcf, merged_vcf,
														geno_map, sample_man,
														option->num_threads);
		if(new_merged_vcf3 != NULL) {
			merged_vcf = new_merged_vcf3;
			continue;
		}
		
		// Impute families whose progenies have been imputed
		auto	*new_merged_vcf4 = impute_vcf_by_progenies(orig_vcf,
													merged_vcf, geno_map,
													sample_man,
													option->num_threads);
		if(new_merged_vcf4 != NULL) {
			merged_vcf = new_merged_vcf4;
			continue;
		}
		break;
	}
	return merged_vcf;
}

VCFSmall *SmallFamily::impute_iolated_samples(
				const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
				SampleManager *sample_man, const STRVEC& samples,
				const Map& gmap, bool modify_genotypes, int num_threads) {
	const STRVEC	references = sample_man->collect_large_family_parents();
	// Split sample to phase for later multithreading
	auto	vcfs_ = VCFIsolated::create(orig_vcf, merged_vcf, samples,
										references, gmap,
										modify_genotypes, num_threads);
	const auto	new_vcfs = VCFIsolated::impute_all(vcfs_, num_threads);
	vector<const VCFSmallBase *>	vcfs(new_vcfs.begin(), new_vcfs.end());
	VCFSmall	*new_vcf = VCFSmall::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(new_vcfs);
	return new_vcf;
}
