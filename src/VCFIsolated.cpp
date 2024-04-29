#include <algorithm>
#include <cassert>
#include "../include/VCFIsolated.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFIsolated ////////////////////

VCFIsolated::VCFIsolated(const vector<STRVEC>& h, const STRVEC& s,
							vector<VCFRecord *> rs, size_t nis,
							const Map& m, bool mg) :
									VCFBase(h, s), VCFImputable(m),
									records(rs), num_imputed_samples(nis),
									modify_genotypes(mg) { }

VCFIsolated::~VCFIsolated() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

vector<Haplotype> VCFIsolated::collect_haplotype_from_refs() const {
	vector<Haplotype>	haps;
	for(size_t i = num_imputed_samples; i < num_samples(); ++i) {
		haps.push_back(clip_haplotype(i, 0));
		haps.push_back(clip_haplotype(i, 1));
	}
	return haps;
}

vector<Haplotype> VCFIsolated::collect_haplotypes_mat(
										size_t sample_index) const {
	return collect_haplotype_from_refs();
}

vector<Haplotype> VCFIsolated::collect_haplotypes_pat(
										size_t sample_index) const {
	return collect_haplotype_from_refs();
}

VCFIsolated *VCFIsolated::divide_by_positions(size_t first, size_t last) const {
	vector<VCFRecord *>	sub_records(records.begin() + first,
									records.begin() + last);
	return new VCFIsolated(get_header(), get_samples(), sub_records,
							num_imputed_samples, get_map(), modify_genotypes);
}

vector<HaplotypePair> VCFIsolated::impute_cM(
									const vector<HaplotypePair>& prev_haps) {
	vector<HaplotypePair>	haps;
	for(size_t i = 0; i < prev_haps.size(); ++i) {
		const auto	prev_hap = prev_haps[i];
		const auto	hap = this->impute_cM_each_sample(prev_hap, i, true,
															modify_genotypes);
		haps.push_back(hap);
	}
	return haps;
}

void VCFIsolated::impute() {
	const Haplotype	h = Haplotype::default_value();
	vector<HaplotypePair>	haps(num_imputed_samples, make_pair(h, h));
	vector<VCFIsolated *>	vcf_cMs = VCFImputable::divide_by_cM(this);
	for(auto p = vcf_cMs.begin(); p != vcf_cMs.end(); ++p) {
		VCFIsolated	*vcf_cM = *p;
		haps = vcf_cM->impute_cM(haps);
		vcf_cM->records.clear();
		delete vcf_cM;
	}
}

VCFSmall *VCFIsolated::extract_isolated_samples() const {
	const STRVEC	isolated_samples(samples.begin(),
									 samples.begin() + num_imputed_samples);
	return extract_samples(isolated_samples);
}

vector<VCFIsolated *> VCFIsolated::create(const VCFSmall *orig_vcf,
											const VCFSmall *imputed_vcf,
											const STRVEC& samples,
											const STRVEC& references,
											const Map& gmap,
											bool modify_genotypes,
											int num_threads) {
	const auto	sample_columns = orig_vcf->extract_columns(samples);
	const auto	ref_columns = imputed_vcf->extract_columns(references);
	const auto	column_table = divide_columns(sample_columns, num_threads);
	vector<VCFIsolated *>	vcfs;
	for(auto p = column_table.begin(); p != column_table.end(); ++p) {
		const auto	cs = *p;
		STRVEC	new_samples;
		for(auto q = cs.begin(); q != cs.end(); ++q)
			new_samples.push_back(orig_vcf->get_samples()[*q-9]);
		new_samples.insert(new_samples.end(),
							references.begin(), references.end());
		const auto	header = orig_vcf->trim_header(new_samples);
		vector<VCFRecord *>	records;
		// Create VCF first and add Record later to have samples in VCF
		VCFIsolated	*vcf = new VCFIsolated(header, new_samples, records,
											cs.size(), gmap, modify_genotypes);
		for(size_t i = 0; i < orig_vcf->size(); ++i) {
			VCFRecord	*record = orig_vcf->get_record(i);
			VCFRecord	*imputed_record = imputed_vcf->get_record(i);
			const STRVEC&	v0 = record->get_v();
			STRVEC	v(v0.begin(), v0.begin() + 9);
			for(auto q = cs.begin(); q != cs.end(); ++q)
				v.push_back(record->get_v()[*q]);
			for(auto q = ref_columns.begin(); q != ref_columns.end(); ++q)
				v.push_back(imputed_record->get_v()[*q]);
			VCFRecord	*new_record = new VCFRecord(v, vcf->get_samples());
			vcf->add_record(new_record);
		}
		vcfs.push_back(vcf);
	}
	return vcfs;
}

vector<vector<size_t>> VCFIsolated::divide_columns(
										const vector<size_t>& cs, int num) {
	vector<vector<size_t>>	css;
	if(cs.size() <= (size_t)num) {
		for(auto p = cs.begin(); p != cs.end(); ++p)
			css.push_back(vector<size_t>(1, *p));
	}
	else {
		css.resize(num);
		for(size_t i = 0; i < cs.size(); ++i)
			css[i%num].push_back(cs[i]);
	}
	return css;
}

vector<VCFSmallBase *> VCFIsolated::impute_all(
									const vector<VCFIsolated *>& vcfs,
									int num_threads) {
	impute_parellel(vcfs, num_threads);
	
	vector<VCFSmallBase *>	new_vcfs;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		new_vcfs.push_back((*p)->extract_isolated_samples());
		delete *p;
	}
	return new_vcfs;
}

void VCFIsolated_impute_in_thread(void *config) {
	const auto	*c = (VCFIsolated::ConfigThread *)config;
	for(size_t i = c->first; i < c->size(); i += c->num_threads) {
		auto	vcf = c->vcfs[i];
		vcf->impute();
	}
}

void VCFIsolated::impute_parellel(const vector<VCFIsolated *>& vcfs,
														int num_threads) {
	const int	T = min((int)vcfs.size(), num_threads);
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(vcfs, (size_t)i, T);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
						(void *(*)(void *))&VCFIsolated_impute_in_thread,
						(void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		VCFIsolated_impute_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
}
