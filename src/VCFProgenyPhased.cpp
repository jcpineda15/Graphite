#include <stack>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <cassert>
#include "../include/common.h"
#include "../include/VCFProgenyPhased.h"
#include "../include/VCFFillable.h"
#include "../include/Map.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/Imputer.h"

using namespace std;


//////////////////// VCFProgenyPhased ////////////////////

VCFProgenyPhased::VCFProgenyPhased(const vector<STRVEC>& h, const STRVEC& s,
									vector<VCFFamilyRecord *> rs,
									const Map& gmap, const VCFSmall *ref) :
				VCFBase(h, s), VCFFamilyBase(), VCFImputable(gmap),
				records(rs), selection(0), ref_vcf(ref) { }

VCFProgenyPhased::~VCFProgenyPhased() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

Haplotype VCFProgenyPhased::clip_ref_haplotype(
								size_t sample_index, int side) const {
	const vector<int>	hap = ref_vcf->clip_raw_haplotype(sample_index, side);
	return Haplotype(hap, sample_index, side);
}

VCFProgenyPhased *VCFProgenyPhased::divide_by_positions(
											size_t first, size_t last) const {
	vector<VCFFamilyRecord *>	rs(records.begin() + first,
								   records.begin() + last);
	vector<VCFRecord *>	ref_rs(ref_vcf->get_records().begin() + first,
							   ref_vcf->get_records().begin() + last);
	const VCFSmall	*ref_vcf_divided = new VCFSmall(ref_vcf->get_header(),
													ref_vcf->get_samples(),
													ref_rs, true);
	return new VCFProgenyPhased(header, samples, rs, get_map(),
														ref_vcf_divided);
}

vector<Haplotype> VCFProgenyPhased::collect_haplotypes_from_phased_progeny(
															int side) const {
	vector<Haplotype>	haps(1, clip_haplotype(2, side));
	return haps;
}

vector<Haplotype> VCFProgenyPhased::collect_haplotype_from_refs() const {
	vector<Haplotype>	haps;
	for(size_t sample_index = 0; sample_index < ref_vcf->num_samples();
														++sample_index) {
		haps.push_back(clip_ref_haplotype(sample_index, 0));
		haps.push_back(clip_ref_haplotype(sample_index, 1));
	}
	return haps;
}

vector<Haplotype> VCFProgenyPhased::collect_haplotypes_mat(
											size_t sample_index) const {
	int side = (this->selection + sample_index) % 2;
	return this->collect_haplotypes_from_phased_progeny(side);
}
vector<Haplotype> VCFProgenyPhased::collect_haplotypes_pat(
											size_t sample_index) const {
	return this->collect_haplotype_from_refs();
}

size_t VCFProgenyPhased::collect_known_parents_indices() const {
	if(this->get_samples()[0] != "0")
		return 0;
	else
		return 1;
}

vector<HaplotypePair> VCFProgenyPhased::impute_core(
									const vector<VCFProgenyPhased *>& vcf_cMs) {
	const size_t	sample_index = vcf_cMs[0]->collect_known_parents_indices();
	
	vector<HaplotypePair>	whole_haplotype_pair = vector<HaplotypePair>();
	const Haplotype	h = Haplotype::default_value();
	HaplotypePair	hap = make_pair(h, h);
	HaplotypePair	*ptr_hap = &hap;
	for(auto p = vcf_cMs.begin(); p != vcf_cMs.end(); ++p) {
		VCFProgenyPhased	*vcf_cM = *p;
		auto	hap_ = vcf_cM->impute_cM_each_sample(*ptr_hap, sample_index,
																false, true);
		whole_haplotype_pair.push_back(hap_);
		ptr_hap = &whole_haplotype_pair.back();
	}
	return whole_haplotype_pair;
}

int VCFProgenyPhased::score_whole(const vector<HaplotypePair>& haps,
							size_t sample_index,
							const vector<VCFProgenyPhased *>& vcf_cMs) const {
	int	score = 0;
	const vector<int>	int_gts = this->get_int_gts(sample_index);
	size_t	first = 0;
	for(size_t i = 0; i < vcf_cMs.size(); ++i) {
		const VCFProgenyPhased *vcf_cM = vcf_cMs[i];
		const HaplotypePair&	hap = haps[i];
		const size_t	last = first + vcf_cM->size();
		const vector<int>	sub_int_gts(int_gts.begin() + first,
										int_gts.begin() + last);
		score += Haplotype::score(hap, sub_int_gts);
		first = last;
	}
	return score;
}

int VCFProgenyPhased::sum_score(const vector<HaplotypePair>& hap_pairs,
								const vector<VCFProgenyPhased *>& vcf_cMs,
								const size_t sample_index) const {
	return this->score_whole(hap_pairs, sample_index, vcf_cMs);
}

void VCFProgenyPhased::impute() {
	const auto	vcf_cMs = VCFImputable::divide_by_cM(this);
	const vector<HaplotypePair>	hap_pairs1 = this->impute_core(vcf_cMs);
	for(auto p = vcf_cMs.begin(); p != vcf_cMs.end(); ++p)
		(*p)->inverse_selection();
	const vector<HaplotypePair>	hap_pairs2 = this->impute_core(vcf_cMs);
	
	const size_t	sample_index = vcf_cMs[0]->collect_known_parents_indices();
	const int	score1 = sum_score(hap_pairs1, vcf_cMs, sample_index);
	const int	score2 = sum_score(hap_pairs2, vcf_cMs, sample_index);
	const auto	hap_pairs = score1 >= score2 ? hap_pairs1 : hap_pairs2;
	for(size_t j = 0; j < vcf_cMs.size(); ++j) {
		vcf_cMs[j]->set_haplotype(hap_pairs[j], sample_index, true);
		vcf_cMs[j]->records.clear();
		delete vcf_cMs[j]->ref_vcf;
		delete vcf_cMs[j];
	}
}

VCFProgenyPhased *VCFProgenyPhased::impute_by_progeny(const VCFSmall *orig_vcf,
													const VCFSmall *imputed_vcf,
													const STRVEC& samples,
													const size_t& ppi,
													const Map& map_,
													const VCFSmall *ref_vcf) {
	VCFProgenyPhased	*vcf = create(orig_vcf, imputed_vcf, samples,
														ppi, map_, ref_vcf);
	vcf->impute();
	return vcf;
}

void VCFProgenyPhased::impute_in_thread(void *config) {
	const auto	*c = (ConfigThread *)config;
	for(size_t i = c->first; i < c->size(); i += c->num_threads) {
		const Family	*family = c->families[i].first;
		const size_t	ppi = c->families[i].second;
		c->results[i] = impute_by_progeny(c->orig_vcf, c->merged_vcf,
												family->get_samples(), ppi,
												c->map, c->ref_vcf);
	}
}

vector<VCFProgenyPhased *> VCFProgenyPhased::impute_all_by_progeny(
						const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
						const vector<PAIR>& families, const Map& gmap,
						const VCFSmall *ref_vcf, int num_threads) {
	vector<VCFProgenyPhased *>	results(families.size());
	
	const int	T = min((int)families.size(), num_threads);
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(orig_vcf, merged_vcf, families,
										gmap, ref_vcf, (size_t)i, T, results);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
						(void *(*)(void *))&impute_in_thread,
						(void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		impute_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	return results;
}

VCFProgenyPhased *VCFProgenyPhased::create(const VCFSmall *orig_vcf,
							const VCFSmall *merged_vcf,
							const STRVEC& family_samples, size_t phased_index,
							const Map& map_, const VCFSmall *ref_vcf) {
	// extract parents and the phased progeny
	const STRVEC	samples { family_samples[0], family_samples[1],
												family_samples[phased_index] };
	VCFFamily	*vcf = VCFFamily::create_by_two_vcfs(merged_vcf,
														orig_vcf, samples);
	auto	vcf2 = new VCFProgenyPhased(vcf->get_header(), vcf->get_samples(),
									vcf->get_family_records(), map_, ref_vcf);
	vcf->clear_records();
	delete vcf;
	return vcf2;
}
