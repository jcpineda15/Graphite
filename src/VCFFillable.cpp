#include <limits>
#include <memory>
#include <cassert>
#include "../include/VCFFillable.h"
#include "../include/VCFHeteroHomo.h"
#include "../include/group.h"

#include "../include/option.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFFillable ////////////////////

VCFFillable::VCFFillable(const std::vector<STRVEC>& h, const STRVEC& s,
								std::vector<VCFFillableRecord *> rs) :
				VCFBase(h, s), VCFSmallBase(), VCFFamilyBase(), records(rs) { }

VCFFillable::~VCFFillable() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

void VCFFillable::phase_in_thread(void *config) {
	auto	*c = (ConfigThreadPhase *)config;
	const auto&	record_sets = c->record_sets;
	const size_t	n = record_sets.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		if(record_sets[i].record->is_fillable_type())
			record_sets[i].impute(true);
	}
}

void VCFFillable::modify(int T) {
	const Groups	*groups = Groups::create(records);
	const auto	record_sets = groups->create_record_sets();
	
	vector<ConfigThreadPhase *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThreadPhase(i, T, record_sets, this);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
			(void *(*)(void *))&phase_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		phase_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	
	delete groups;
	
	for(size_t i = 0U; i < records.size(); ++i) {
		auto	*record = records[i];
		if(record->is_mat_type())
			this->impute_NA_mat(i);
		else if(record->is_pat_type())
			this->impute_NA_pat(i);
		else if(record->is_fillable_type())
			this->impute_others(i);
	}
	
	for(auto p = records.begin(); p != records.end(); ++p)
		(*p)->fill_PGT();
}

void VCFFillable::phase_hetero_hetero() {
	// Group records with type 'IMPUTABLE', 'MAT', 'PAT' and 'FIXED'
	const Groups	*groups = Groups::create(records);
	const auto	record_sets = groups->create_record_sets();
	for(auto p = record_sets.begin(); p != record_sets.end(); ++p) {
		p->impute(false);
	}
	
	for(size_t i = 0U; i < records.size(); ++i) {
		auto	*record = records[i];
		if(record->is_mat_type())
			this->impute_NA_mat(i);
		else if(record->is_pat_type())
			this->impute_NA_pat(i);
		else if(record->is_fillable_type())
			this->impute_others(i);
	}
	delete groups;
}

VCFFillable *VCFFillable::create_from_header() const {
	vector<VCFFillableRecord *>	rs;
	VCFFillable	*vcf = new VCFFillable(header, samples, rs);
	copy_chrs(vcf);
	return vcf;
}

template<typename Iter>
VCFFillableRecord *VCFFillable::find_neighbor_same_type_record(
							size_t i, size_t c, Iter first, Iter last) const {
	const FillType	type = records[i]->get_type();
	const string&	chromosome = records[i]->chrom();
	for(auto p = first; p != last; ++p) {
		auto	*record = *p;
		if(record->chrom() != chromosome)
			return NULL;
		else if(record->get_type() == type && record->get_GT(c-9) != "./.")
			return record;
	}
	return NULL;
}

VCFFillableRecord *VCFFillable::find_prev_same_type_record(
													size_t i, size_t c) const {
	if(i == 0U)
		return NULL;
	
	return find_neighbor_same_type_record(i, c, records.rend() - i + 1,
																records.rend());
}

VCFFillableRecord *VCFFillable::find_next_same_type_record(
													size_t i, size_t c) const {
	if(i == records.size() - 1)
		return NULL;
	
	const auto	first = records.begin() + i + 1;
	const auto	last = records.end();
	return find_neighbor_same_type_record(i, c, first, last);
}

const RecordSet *VCFFillable::create_recordset(
										size_t i, size_t c, bool is_mat) const {
	auto	*record = records[i];
	auto	*prev_record = find_prev_same_type_record(i, c);
	auto	*next_record = find_next_same_type_record(i, c);
	if(is_mat)
		return new RecordSet(record, prev_record, next_record, NULL, NULL);
	else
		return new RecordSet(record, NULL, NULL, prev_record, next_record);
}

void VCFFillable::impute_NA_mat_each(size_t i, size_t c) {
	const RecordSet	*rs = this->create_recordset(i, c, true);
	rs->impute_NA_mat_each(c - 9);
	delete rs;
}

// impute N/A genotypes in families with heterozygous mat and homozygous pat
void VCFFillable::impute_NA_mat(size_t i) {
	auto	*record = this->records[i];
	for(size_t c = 11U; c != samples.size() + 9; ++c) {
		if(record->get_GT(c-9) == "./.")
			impute_NA_mat_each(i, c);
	}
}

void VCFFillable::impute_NA_pat_each(size_t i, size_t c) {
	const RecordSet	*rs = this->create_recordset(i, c, false);
	rs->impute_NA_pat_each(c - 9);
	delete rs;
}

void VCFFillable::impute_NA_pat(size_t i) {
	auto	*record = this->records[i];
	for(size_t c = 11U; c != samples.size() + 9; ++c) {
		if(record->get_GT(c-9) == "./.")
			impute_NA_pat_each(i, c);
	}
}

pair<int, int> VCFFillable::find_prev_mat_from(int i, int c) const {
	for(int k = i - 1; k >= 0; --k) {
		const int	from1 = this->records[k]->mat_from(c);
		if(from1 != 0)
			return pair<int, int>(k, from1);
	}
	return pair<int, int>(-1, 0);
}

pair<int, int> VCFFillable::find_next_mat_from(int i, int c) const {
	for(int k = i + 1; k < (int)this->size(); ++k) {
		const int	from2 = this->records[k]->mat_from(c);
		if(from2 != 0)
			return pair<int, int>(k, from2);
	}
	return pair<int, int>(-1, 0);
}

pair<int, int> VCFFillable::find_prev_pat_from(int i, int c) const {
	for(int k = i - 1; k >= 0; --k) {
		const int	from1 = this->records[k]->pat_from(c);
		if(from1 != 0)
			return pair<int, int>(k, from1);
	}
	return pair<int, int>(-1, 0);
}

pair<int, int> VCFFillable::find_next_pat_from(int i, int c) const {
	for(int k = i + 1; k < (int)this->size(); ++k) {
		const int	from2 = this->records[k]->pat_from(c);
		if(from2 != 0)
			return pair<int, int>(k, from2);
	}
	return pair<int, int>(-1, 0);
}

int VCFFillable::select_from(const pair<int, int>& f1,
							 const pair<int, int>& f2, int i) const {
	const int	i1 = f1.first;
	const int	from1 = f1.second;
	const int	i2 = f2.first;
	const int	from2 = f2.second;
	if(from1 == 0 && from2 == 0) {
		// Decide which Haplotype comes from when there is no before or after
		const auto	*r0 = this->records[i];
		return r0->pos() % 2 + 1;
	}
	if(from1 == from2)
		return from1;
	else if(from2 == 0)
		return from1;
	else if(from1 == 0)
		return from2;
	else {
		// Determine which Haplotype is coming from by physical distance
		const auto	*r0 = this->records[i];
		const auto	*r1 = this->records[i1];
		const auto	*r2 = this->records[i2];
		if(r0->pos() * 2 <= r1->pos() + r2->pos())
			return from1;
		else
			return from2;
	}
}

int VCFFillable::find_mat_from(int i, int c) const {
	return this->select_from(this->find_prev_mat_from(i, c),
							 this->find_next_mat_from(i, c), i);
}

int VCFFillable::find_pat_from(int i, int c) const {
	return this->select_from(this->find_prev_pat_from(i, c),
							 this->find_next_pat_from(i, c), i);
}

void VCFFillable::impute_others(int i) {
	auto	*record = this->records[i];
	const bool	mat_homo = record->is_homo(0);
	const bool	pat_homo = record->is_homo(1);
	for(size_t c = 11; c != record->get_v().size(); ++c) {
		if(record->get_v()[c].c_str()[1] != '/' &&
								record->get_int_gt(c-9) != -1)
			continue;
		const int	mat_from = mat_homo ? 1 : this->find_mat_from(i, c);
		const int	pat_from = pat_homo ? 1 : this->find_pat_from(i, c);
		record->set_GT(c-9, record->gt_from_parent(mat_from, pat_from));
	}
}

VCFFillable *VCFFillable::fill(const vector<VCFHeteroHomo *>& vcfs,
				const vector<VCFImpFamilyRecord *>& records, int num_threads) {
	vector<VCFFillableRecord *>	merged_records
							 = VCFFillable::merge_records(vcfs, records, true);
	VCFFillable	*vcf = new VCFFillable(vcfs.front()->get_header(),
								vcfs.front()->get_samples(), merged_records);
	vcf->modify(num_threads);
	return vcf;
}

vector<VCFFillableRecord *> VCFFillable::merge_records(
									const vector<VCFHeteroHomo *>& vcfs,
									const vector<VCFImpFamilyRecord *>& records,
									bool all_out) {
	vector<VCFFillableRecord *>	all_records;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		auto	*vcf = *p;
		auto&	hh_records = vcf->get_records();
		for(auto q = hh_records.begin(); q != hh_records.end(); ++q)
			all_records.push_back(VCFFillableRecord::convert(*q));
	}
	
	for(auto p = records.begin(); p != records.end(); ++p) {
		if(all_out || (*p)->is_fillable())
			all_records.push_back(VCFFillableRecord::convert(*p));
	}
	std::sort(all_records.begin(), all_records.end(),
				[](const VCFFillableRecord *r1, const VCFFillableRecord *r2) {
					return r1->get_index() < r2->get_index(); });
	return all_records;
}

vector<vector<VCFFillableRecord *>> VCFFillable::collect_records(
										const vector<VCFFillable *>& vcfs) {
	vector<vector<VCFFillableRecord *>>	rss;
	for(size_t i = 0; i < vcfs.front()->size(); ++i) {
		vector<VCFFillableRecord *>	rs;
		for(auto p = vcfs.begin(); p != vcfs.end(); ++p)
			rs.push_back((*p)->get_fillable_record(i));
		rss.push_back(rs);
	}
	return rss;
}

pair<STRVEC, vector<vector<pair<int, int>>>>
VCFFillable::integrate_samples(const vector<STRVEC>& sss,
									const STRVEC& orig_samples) {
	// { sample: (index of family, index of inner family) }
	map<string, vector<pair<int, int>>>	dic;
	for(size_t i = 0; i < sss.size(); ++i) {
		const vector<string>&	samples = sss[i];
		for(size_t j = 0; j < samples.size(); ++j)
			dic[samples[j]].push_back(make_pair(i, j));
	}
	
	STRVEC	new_samples;
	vector<vector<pair<int, int>>>	pos_samples;
	for(auto p = orig_samples.begin(); p != orig_samples.end(); ++p) {
		const string&	sample = *p;
		if(dic.find(sample) != dic.end()) {
			new_samples.push_back(sample);
			pos_samples.push_back(dic[sample]);
		}
	}
	
	return make_pair(new_samples, pos_samples);
}

// Integrate the VCF so that duplicate samples are one
VCFSmall *VCFFillable::integrate(const VCFFillable *vcf,
								const vector<vector<VCFFillableRecord *>>& rss,
								const STRVEC& orig_samples) {
	vector<STRVEC>	sss;
	for(auto p = rss.front().begin(); p != rss.front().end(); ++p) {
		const STRVEC&	samples = (*p)->get_samples();
		sss.push_back(samples);
	}
	
	const auto	q = integrate_samples(sss, orig_samples);
	const STRVEC&	new_samples = q.first;
	// What number of each Family is the sample at?
	const auto&		pos_samples = q.second;
	
	const vector<STRVEC>	header = vcf->trim_header(new_samples);
	VCFSmall	*new_vcf = new VCFSmall(header, new_samples,
											vector<VCFRecord *>());
	const STRVEC&	samples = new_vcf->get_samples();
	vector<VCFRecord *>	records;
	for(auto p = rss.begin(); p != rss.end(); ++p) {
		auto	*r = VCFFillableRecord::integrate(*p, samples, pos_samples);
		records.push_back(r);
	}
	new_vcf->add_records(records);
	return new_vcf;
}

VCFSmall *VCFFillable::merge(const vector<VCFFillable *>& vcfs,
											const STRVEC& orig_samples) {
	const auto	rss = collect_records(vcfs);
	return integrate(vcfs.front(), rss, orig_samples);
}

void VCFFillable::fill_in_thread(void *config) {
	const auto	*c = (ConfigFillThread *)config;
	const size_t	n = c->size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		auto	vcfs = c->items[i].first;
		auto	records = c->items[i].second;
		auto	result = VCFFillable::fill(vcfs, records, c->num_threads);
		c->filled_vcfs[i] = result;
	}
}

vector<VCFFillable *> VCFFillable::fill_parellel(vector<Item>& items,
														int num_threads) {
	vector<VCFFillable *>	results(items.size());
	
	const int	T = min((int)items.size(), num_threads);
	vector<ConfigFillThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigFillThread(items, true, i, T, results);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
					(void *(*)(void *))&fill_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		fill_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	
	return results;
}

void VCFFillable::delete_items(const vector<Item>& items) {
	for(auto q = items.begin(); q != items.end(); ++q) {
		for(auto r = q->first.begin(); r != q->first.end(); ++r)
			delete *r;
		for(auto r = q->second.begin(); r != q->second.end(); ++r)
			delete *r;
	}
}

vector<VCFFillable *> VCFFillable::fill_all(
							map<Parents, vector<VCFHeteroHomo *>>& imputed_vcfs,
							ImpRecords& other_records, int num_threads) {
	vector<Item>	items;
	for(auto q = imputed_vcfs.begin(); q != imputed_vcfs.end(); ++q) {
		const Parents&	parents = q->first;
		vector<VCFHeteroHomo *>&	vcfs = q->second;
		items.push_back(make_pair(vcfs, other_records[parents]));
	}
	vector<VCFFillable *>	filled_vcfs = fill_parellel(items, num_threads);
	delete_items(items);
	return filled_vcfs;
}
