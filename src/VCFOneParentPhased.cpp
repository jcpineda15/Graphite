#include <stack>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <cassert>
#include "../include/common.h"
#include "../include/VCFOneParentPhased.h"
#include "../include/VCFFillable.h"
#include "../include/Map.h"
#include "../include/Pedigree.h"
#include "../include/Imputer.h"

using namespace std;


//////////////////// VCFOneParentPhased ////////////////////

VCFOneParentPhased::VCFOneParentPhased(const vector<STRVEC>& h, const STRVEC& s,
									   vector<VCFFamilyRecord *> rs, bool mat_p,
									   const Map& m, const VCFSmall *ref) :
						VCFBase(h, s), VCFFamilyBase(), VCFImputable(m),
						records(rs), is_mat_phased(mat_p), ref_vcf(ref) { }

VCFOneParentPhased::~VCFOneParentPhased() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

bool VCFOneParentPhased::is_mat_hetero() const {
	return records.front()->is_hetero(0);
}

Haplotype VCFOneParentPhased::clip_haplotype(
									size_t sample_id, int i) const {
	const auto	hap = clip_raw_haplotype(sample_id, i);
	return Haplotype(hap, sample_id, i);
}

Haplotype VCFOneParentPhased::clip_ref_haplotype(
									size_t sample_id, int i) const {
	const auto	hap = ref_vcf->clip_raw_haplotype(sample_id, i);
	return Haplotype(hap, sample_id, i);
}

VCFOneParentPhased *VCFOneParentPhased::divide_by_positions(
											size_t first, size_t last) const {
	vector<VCFFamilyRecord *>	rs(records.begin() + first,
								   records.begin() + last);
	vector<VCFRecord *>	ref_rs(ref_vcf->get_records().begin() + first,
							   ref_vcf->get_records().begin() + last);
	const VCFSmall	*ref_vcf_divided = new VCFSmall(ref_vcf->get_header(),
													ref_vcf->get_samples(),
													ref_rs, true);
	return new VCFOneParentPhased(header, samples, rs, is_mat_phased, get_map(),
															ref_vcf_divided);
}

vector<Haplotype> VCFOneParentPhased::collect_haplotypes_from_parents() const {
	const size_t	i = this->is_mat_phased ? 0 : 1;
	const vector<Haplotype>	haps { clip_haplotype(i, 0), clip_haplotype(i, 1) };
	return haps;
}

vector<Haplotype> VCFOneParentPhased::collect_haplotype_from_refs() const {
	vector<Haplotype>	haps;
	for(size_t i = 0; i < ref_vcf->num_samples(); ++i) {
		haps.push_back(clip_ref_haplotype(i, 0));
		haps.push_back(clip_ref_haplotype(i, 1));
	}
	return haps;
}

vector<Haplotype> VCFOneParentPhased::collect_haplotypes_mat(
												size_t sample_index) const {
	if(this->is_mat_phased)
		return collect_haplotypes_from_parents();
	else
		return collect_haplotype_from_refs();
}

vector<Haplotype> VCFOneParentPhased::collect_haplotypes_pat(
												size_t sample_index) const {
	if(this->is_mat_phased)
		return collect_haplotype_from_refs();
	else
		return collect_haplotypes_from_parents();
}

vector<HaplotypePair> VCFOneParentPhased::impute_cM(
									const vector<HaplotypePair>& prev_haps) {
	vector<HaplotypePair>	haps;
	for(size_t i = 0; i < prev_haps.size(); ++i) {
		const auto	prev_hap = prev_haps[i];
		const auto	hap = this->impute_cM_each_sample(prev_hap, i+2,
															true, true);
		haps.push_back(hap);
	}
	return haps;
}

void VCFOneParentPhased::impute() {
	const Haplotype	h = Haplotype::default_value();
	vector<HaplotypePair>	haps(num_progenies(), make_pair(h, h));
	const auto	vcf_cMs = VCFImputable::divide_by_cM(this);
	for(auto p = vcf_cMs.begin(); p != vcf_cMs.end(); ++p) {
		VCFOneParentPhased	*vcf_cM = *p;
		haps = vcf_cM->impute_cM(haps);
		vcf_cM->records.clear();
		delete vcf_cM->ref_vcf;
		delete vcf_cM;
	}
}

VCFOneParentPhased *VCFOneParentPhased::create(
						const STRVEC& samples, bool is_mat_phased,
						const VCFSmall *merged_vcf, const VCFSmall *orig_vcf,
						const Map& gmap, const VCFSmall *ref_vcf) {
	VCFFamily	*vcf = VCFFamily::create_by_two_vcfs(merged_vcf,
														orig_vcf, samples);
	vector<VCFFamilyRecord *>	records = vcf->get_family_records();
	auto	new_vcf = new VCFOneParentPhased(vcf->get_header(),
												vcf->get_samples(), records,
												is_mat_phased, gmap, ref_vcf);
	vcf->clear_records();
	delete vcf;
	return new_vcf;
}

void VCFOneParentPhased::impute_in_thread(void *config) {
	const auto	*c = (VCFOneParentPhased::ConfigThread *)config;
	for(size_t i = c->first; i < c->size(); i += c->num_threads) {
		VCFOneParentPhased	*vcf = c->vcfs[i];
		vcf->impute();
	}
}

void VCFOneParentPhased::impute_in_parallel(
					const vector<VCFOneParentPhased *>& vcfs, int num_threads) {
	const int	T = min((int)vcfs.size(), num_threads);
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(vcfs, (size_t)i, T);
	
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
}

STRVEC VCFOneParentPhased::collect_samples(
								const vector<VCFOneParentPhased *>& vcfs) {
	set<string>	set_samples;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const auto	*vcf = *p;
		if(vcf->is_mat_phased)
			set_samples.insert(vcf->pat());
		else
			set_samples.insert(vcf->mat());
		for(size_t i = 2; i < vcf->num_samples(); ++i)
			set_samples.insert(vcf->samples[i]);
	}
	return STRVEC(set_samples.begin(), set_samples.end());
}

// make a record corresponding to samples
VCFSmall *VCFOneParentPhased::merge(const vector<VCFOneParentPhased *>& vcfs) {
	const STRVEC	samples = collect_samples(vcfs);
	
	// Where sample is located in which VCF
	map<string, pair<const VCFOneParentPhased *, size_t>>	dic_pos;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const VCFOneParentPhased	*vcf = *p;
		for(size_t k = 0; k < vcf->num_samples(); ++k)
			dic_pos[vcf->get_samples()[k]] = make_pair(vcf, k + 9);
	}
	
	// create records from the above dictionary
	vector<VCFRecord *>	records;
	for(size_t i = 0; i < vcfs.front()->size(); ++i) {
		const STRVEC&	v1 = vcfs.front()->get_record(i)->get_v();
		STRVEC	v(v1.begin(), v1.begin() + 9);
		for(auto p = samples.begin(); p != samples.end(); ++p) {
			auto	q = dic_pos[*p];
			v.push_back(q.first->get_record(i)->get_v()[q.second]);
		}
		records.push_back(new VCFRecord(v, samples));
	}
	vector<STRVEC>	header = vcfs.front()->trim_header(samples);
	return new VCFSmall(header, samples, records);
}

VCFSmall *VCFOneParentPhased::impute_all(
					const vector<VCFOneParentPhased *>& vcfs, int num_threads) {
	impute_in_parallel(vcfs, num_threads);
	return merge(vcfs);
}
