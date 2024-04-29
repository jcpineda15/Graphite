#include <sstream>
#include <cassert>

#include "../include/VCFFillableRecord.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFFillableRecord ////////////////////

#define P(a,b) make_pair(a,b)

vector<pair<int, int>> VCFFillableRecord::possible_phasings() const {
	switch(this->comb) {
		case ParentComb::P00x00: return { P(0, 0) };
		case ParentComb::P00x01: return { P(0, 1), P(0, 2), P(1, 0), P(2, 0),
										  P(1, 1), P(1, 2), P(2, 1), P(2, 2) };
		case ParentComb::P01x01: return { P(1, 1), P(1, 2), P(2, 1), P(2, 2),
										  P(0, 1), P(0, 2), P(1, 0), P(2, 0),
										  P(1, 3), P(2, 3), P(3, 1), P(3, 2) };
		case ParentComb::P00x11: return { P(0, 3), P(3, 0) };
		case ParentComb::P01x11: return { P(1, 3), P(2, 3), P(3, 1), P(3, 2),
										  P(1, 1), P(1, 2), P(2, 1), P(2, 2) };
		case ParentComb::P11x11: return { P(3, 3) };
		default: break;
	}
	
	// any pair is possible
	vector<pair<int, int>>	pairs;
	for(int i = 0; i < 16; ++i)
		pairs.push_back(P(i >> 2, i & 3));
	return pairs;
}

STRVEC VCFFillableRecord::prog_gts() const {
	STRVEC	vec(v.size() - 11);
	std::copy(v.begin() + 11, v.end(), vec.begin());
	return vec;
}

VCFFillableRecord *VCFFillableRecord::copy() const {
	return new VCFFillableRecord(v, samples, index, type, comb);
}

tuple<int,int,int> VCFFillableRecord::count_gts() const {
	const auto	gts = get_progeny_int_gts();
	int	counter[3] = { 0, 0, 0 };
	for(auto p = gts.begin() + 2; p != gts.end(); ++p) {
		if(*p != -1)
			counter[*p] += 1;
	}
	return make_tuple(counter[0], counter[1], counter[2]);
}

string VCFFillableRecord::gt_from_parent(int mat_from, int pat_from) const {
	stringstream	ss;
	ss << v[9].c_str()[mat_from*2-2] << '|' << v[10].c_str()[pat_from*2-2];
	return ss.str();
}

string VCFFillableRecord::gt_from_mat(int mat_from, int c) const {
	const string	gt = this->get_GT(c - 9);
	const int		int_gt = this->get_int_gt(c - 9);
	const char		mat_gt = this->v[9].c_str()[mat_from*2-2];
	const string	pat_GT = this->get_GT(1);
	if(int_gt == 0) {
		if(mat_gt == '0')
			return "0|0";
		else if(pat_GT.find('0') != string::npos)
			return "1|0";
		else
			return "1|1";
	}
	else if(int_gt == 2) {
		if(mat_gt == '1')
			return "1|1";
		else if(pat_GT.find('1') != string::npos)
			return "0|1";
		else
			return "0|0";
	}
	else if(int_gt == 1) {
		if(mat_gt == '0')
			return "0|1";
		else
			return "1|0";
	}
	else {
		return ".|.";
	}
}

string VCFFillableRecord::gt_from_pat(int pat_from, int c) const {
	const string	gt = this->get_GT(c - 9);
	const int		int_gt = this->get_int_gt(c - 9);
	const char		pat_gt = this->v[10].c_str()[pat_from*2-2];
	const string	mat_GT = this->get_GT(0);
	if(int_gt == 0) {
		if(pat_gt == '0')
			return "0|0";
		else if(mat_GT.find('0') != string::npos)
			return "0|1";
		else
			return "1|1";
	}
	else if(int_gt == 2) {
		if(pat_gt == '1')
			return "1|1";
		else if(mat_GT.find('1') != string::npos)
			return "1|0";
		else
			return "0|0";
	}
	else if(int_gt == 1) {
		if(pat_gt == '0')
			return "1|0";
		else
			return "1|0";
	}
	else {
		return ".|.";
	}
}

vector<vector<double>> VCFFillableRecord::make_probability_table() const {
	// Probability for all parent Genotype combinations
	// (0/0, 0/0), (0/0, 0/1), (0/0, 1/1), (0/1, 0/1), (0/1, 1/1), (1/1, 1/1)
	const vector<double>	ps = {1.0, 0.5, 0.0};
	vector<vector<double>>	pss_;
	for(int i = 0; i < 3; ++i) {
		const double	mat = ps[i];
		for(int j = i; j < 3; ++j) {
			const double	pat = ps[j];
			const vector<double>	v = {mat*pat, mat+pat-2*mat*pat,
													(1.0-mat)*(1.0-pat)};
			pss_.push_back(v);
		}
	}
	
	vector<vector<double>>	pss(6, vector<double>(3));
	const double	p_miss = 0.01;
	for(int i = 0; i < 6; ++i) {
		for(int j = 0; j < 3; ++j)
			pss[i][j] = (pss_[i][j] + p_miss) / (1.0 + 3 * p_miss);
	}
	return pss;
}

void VCFFillableRecord::phase() {
	// premise that both parents are homozygous
	const string	gt_mat = this->v[9].substr(0, 1);
	const string	gt_pat = this->v[10].substr(0, 1);
	this->set_mat_GT(gt_mat + "|" + gt_mat);
	this->set_pat_GT(gt_pat + "|" + gt_pat);
	for(size_t i = 2; i < this->num_samples(); ++i)
		this->set_GT(i, gt_mat + "|" + gt_pat);
}

int VCFFillableRecord::mat_from(int c) const {
	if(this->v[c].c_str()[0] == '.')
		return 0;
	else if(!this->is_hetero(0))
		return 0;
	else if(this->v[c].c_str()[0] == this->v[9].c_str()[0])
		return 1;
	else
		return 2;
}

int VCFFillableRecord::pat_from(int c) const {
	if(this->v[c].c_str()[0] == '.')
		return 0;
	else if(!this->is_hetero(1))
		return 0;
	else if(this->v[c].c_str()[2] == this->v[10].c_str()[0])
		return 1;
	else
		return 2;
}

int VCFFillableRecord::find_geno_type(const string& type) const {
	const auto	vec = Common::split(this->v[8], ':');
	for(auto p = vec.begin(); p != vec.end(); ++p) {
		if(*p == type)
			return p - vec.begin();
	}
	return -1;
}

void VCFFillableRecord::fill_PGT() {
	const int	i_GT = this->find_geno_type("GT");
	assert(i_GT != -1);
	const int	i_PGT = this->find_geno_type("PGT");
	if(i_PGT == -1)
		return;
	
	for(size_t j = 9U; j < this->v.size(); ++j) {
		auto	vec = Common::split(this->v[j], ':');
		vec[i_PGT] = v[i_GT];
		this->v[j] = Common::join(vec, ':');
	}
}

string VCFFillableRecord::inverse_gt(const string& gt, bool inv) const {
	if(inv)
		return string(gt == "0" ? "1" : "0");
	else
		return gt;
}

string VCFFillableRecord::inverse_prog_gt(const string& gt,
											bool inv_mat, bool inv_pat) const {
	return inverse_gt(gt.substr(0, 1), inv_mat) + "|" +
			inverse_gt(gt.substr(2, 1), inv_pat) + gt.substr(3);
}

STRVEC VCFFillableRecord::inverse_prog_gts(const STRVEC& prog_gts,
											bool inv_mat, bool inv_pat) const {
	STRVEC	inv_prog_gts;
	for(auto p = prog_gts.begin(); p != prog_gts.end(); ++p) {
		inv_prog_gts.push_back(inverse_prog_gt(*p, inv_mat, inv_pat));
	}
	return inv_prog_gts;
}

void VCFFillableRecord::inverse_parents_gts(bool inv_mat, bool inv_pat) {
	// both must be hetero
	if(inv_mat)
		set_GT(0, v[9].substr(2, 1) + "|" + v[9].substr(0, 1));
	if(inv_pat)
		set_GT(1, v[10].substr(2, 1) + "|" + v[10].substr(0, 1));
}

bool VCFFillableRecord::is_same_gts(const string& gt1,
									const string& gt2) const {
	if(gt2 == "0/1")
		return gt1 == "0|1" || gt1 == "1|0";
	else if(gt2 == "0/0")
		return gt1 == "0|0";
	else if(gt2 == "1/1")
		return gt1 == "1|1";
	else
		return false;
}

bool VCFFillableRecord::is_near_prog_gts(const STRVEC& gts) const {
	int	num = 0;
	int	dist = 0;
	for(size_t i = 0U; i < gts.size(); ++i) {
		if(v[i+11] != "0/1")
			num += 1;
		if(!is_same_gts(gts[i], v[i+11]))
			dist += 1;
	}
	return dist < max(1, num / 2);
}

void VCFFillableRecord::modify_gts(const STRVEC& new_prog_gts) {
	for(int i = 0; i < 4; ++i) {
		const bool	inv_mat = (i >> 1) == 1;
		const bool	inv_pat = (i & 1) == 1;
		STRVEC	inv_prog_gts = inverse_prog_gts(new_prog_gts, inv_mat, inv_pat);
		if(is_near_prog_gts(inv_prog_gts)) {
			inverse_parents_gts(inv_mat, inv_pat);
			std::copy(inv_prog_gts.begin(), inv_prog_gts.end(), v.begin() + 11);
			return;
		}
	}
	std::copy(new_prog_gts.begin(), new_prog_gts.end(), v.begin() + 11);
}

void VCFFillableRecord::modify_parents_type() {
	if(!is_00x11() &&
			((get_GT(0) == "0|0" && get_GT(1) == "1|1") ||
			 (get_GT(0) == "1|1" && get_GT(1) == "0|0")))
		comb = ParentComb::P00x11;
}

int VCFFillableRecord::from_which_chrom(size_t i, bool is_mat) const {
	int j = is_mat ? 0 : 1;
	const string&	parent_gt = this->get_gt(j);
	return parent_gt.c_str()[0] == this->get_gt(i).c_str()[j*2] ? 1 : 2;
}

void VCFFillableRecord::set(const STRVEC& new_v, FillType new_type) {
	v = new_v;
	type = new_type;
	VCFFamilyRecord::set(new_v);
}

VCFRecord *VCFFillableRecord::integrate_records(
								const vector<VCFFillableRecord *>& records) {
	for(auto p = records.begin(); p != records.end(); ++p) {
		if((*p)->is_unable())
			return NULL;
	}
	
	STRVEC	v(records.front()->v.begin(), records.front()->v.begin() + 9);
	for(auto p = records.begin(); p != records.end(); ++p) {
		v.insert(v.end(), (*p)->v.begin() + 9, (*p)->v.end());
	}
	return new VCFRecord(v, records.front()->get_samples());
}

int VCFFillableRecord::hash(int d) const {
	int	hash_ = 0;
	const string&	str_pos = this->v[1];
	for(auto p = str_pos.begin(); p != str_pos.end(); ++p) {
		hash_ = (hash_ + (*p - '0')) % d;
	}
	return hash_;
}

string VCFFillableRecord::decide_by_majority(const vector<string>& GTs) const {
	// ex) ['0/0', '0/1', '0/0', '1/1'] -> { '0/0': 2, '0/1': 1, '1/1': 1 }
	map<string, int>	counter;
	for(auto p = GTs.begin(); p != GTs.end(); ++p)
		counter[*p] += 1;
	
	// ex) { '0/0': 2, '0/1': 1, '1/1': 1 } -> { 2: ['0/0'], 1: ['0/1', '1/1'] }
	map<int, vector<string>>	dic_num;
	for(auto p = counter.begin(); p != counter.end(); ++p)
		dic_num[p->second].push_back(p->first);
	
	pair<int, vector<string>>	max_pair;
	for(auto p = dic_num.begin(); p != dic_num.end(); ++p) {
		if(p->first > max_pair.first)
			max_pair = *p;
	}
	
	vector<string>	max_GTs = max_pair.second;
	std::sort(max_GTs.begin(), max_GTs.end());	// to match the Python version
	if(max_GTs.size() == 1)
		return max_GTs.front();
	else	// Use a random number if there is more than one candidates
		return max_GTs[this->hash((int)max_GTs.size())];
}

void VCFFillableRecord::swap_parents(int i, const string& GT) {
	if(GT == this->v[i+9].substr(0, 3))
		return;
	
	if(GT == "0|0" || GT == "1|1") {
		const bool	is_mat_00 = (i == 0) ^ (GT == "1|1");
		this->set_GT(0, is_mat_00 ? "0|0" : "1|1");
		this->set_GT(1, is_mat_00 ? "1|1" : "0|0");
		
		// swap genotyeps of progenies
		const string	prog_GT = is_mat_00 ? "0|1" : "1|0";
		for(size_t i = 2; i < this->samples.size(); ++i)
			this->set_GT(i, prog_GT);
	}
}

// Determine a single Genotype
// when the Genotype of the same sample differs from family to family
string VCFFillableRecord::decide_duplicated_Genotype(
									const vector<VCFFillableRecord *>& records,
									const vector<pair<int, int>>& positions) {
	STRVEC	GTs;
	for(auto p = positions.begin(); p != positions.end(); ++p)
		GTs.push_back(records[p->first]->get_GT(p->second));
	
	// if all genotype are ./., not change
	if(Genotype::is_all_NA(GTs))
		return "./.";
	
	// if the sample is a progeny of a family, use the genotype
	for(size_t i = 0; i < GTs.size(); ++i) {
		const int	j = positions[i].second;	// jth(0-based) sample of family
		const string&	GT = GTs[i];
		if(j >= 2 && GT != "./.")	// j >= 2 means progeny
			return GT;
	}
	
	STRVEC	GTs_less_NA;
	for(auto p = GTs.begin(); p != GTs.end(); ++p) {
		if(*p != "./.")
			GTs_less_NA.push_back(*p);
	}
	
	// If all Genotypes except ./. are the same Genotype, use that Genotype
	if(Common::is_all_same(GTs_less_NA))
		return GTs_less_NA.front();
	
	STRVEC	GTs_less_00x11;
	for(size_t i = 0; i < GTs.size(); ++i) {
		const auto	*record = records[positions[i].first];
		const string&	GT = GTs[i];
		if(GT != "./." && record->comb != ParentComb::P00x11)
			GTs_less_00x11.push_back(GT);
	}
	
	if(GTs_less_00x11.empty())	// only 0/0 x 1/1
		return records.front()->decide_by_majority(GTs_less_NA);
	else if(Common::is_all_same(GTs_less_00x11))
		return GTs_less_00x11.front();
	else
		return records.front()->decide_by_majority(GTs_less_00x11);
}

VCFFillableRecord *VCFFillableRecord::convert(
									const VCFImpFamilyRecord *record) {
	const FillType	type = record->get_fill_type();
	return new VCFFillableRecord(record->get_v(), record->get_samples(),
								 record->get_index(), type, record->get_comb());
}

VCFRecord *VCFFillableRecord::merge(const vector<VCFFillableRecord *>& records,
														const STRVEC& samples) {
	STRVEC	v = records.front()->v;
	for(auto p = records.begin() + 1; p != records.end(); ++p)
		v.insert(v.end(), (*p)->v.begin(), (*p)->v.end());
	return new VCFRecord(v, samples);
}

void VCFFillableRecord::integrate_each_sample(
									const vector<VCFFillableRecord *>& records,
									const vector<pair<int, int>>& positions) {
	const string	GT = decide_duplicated_Genotype(records, positions);
	for(auto p = positions.begin(); p != positions.end(); ++p) {
		VCFFillableRecord	*record = records[p->first];
		if(record->is_00x11() && p->second <= 1)
			record->swap_parents(p->second, GT);
	}
}

bool VCFFillableRecord::is_all_same_GT(
						const vector<VCFFillableRecord *>& records,
						const vector<pair<int, int>>& pos_samples) {
	if(pos_samples.size() == 1)
		return true;
	
	const auto	pos_sample = pos_samples.front();
	const string	GT0 = records[pos_sample.first]->get_GT(pos_sample.second);
	for(auto p = pos_samples.begin() + 1; p != pos_samples.end(); ++p) {
		const string	GT = records[p->first]->get_GT(p->second);
		if(GT != GT0)
			return false;
	}
	return true;
}

VCFRecord *VCFFillableRecord::integrate(
							const vector<VCFFillableRecord *>& records,
							const vector<string>& samples,
							const vector<vector<pair<int, int>>>& pos_samples) {
	const VCFFillableRecord	*record = records.front();
	// Exchangeable if one parent Genotype is 0|0 and the other is 1|1
	for(auto p = pos_samples.begin(); p != pos_samples.end(); ++p) {
		if(!is_all_same_GT(records, *p)) {
			integrate_each_sample(records, *p);
		}
	}
	
	vector<string>	v(record->v.begin(), record->v.begin() + 9);
	for(auto p = pos_samples.begin(); p != pos_samples.end(); ++p) {
		const pair<int, int>&	pos = p->front();
		v.push_back(records[pos.first]->get_gt(pos.second));
	}
	return new VCFRecord(v, samples);
}

// after imputed
int VCFFillableRecord::from_which_chrom(const VCFFillableRecord *record,
													size_t i, bool is_mat) {
	if(record == NULL)
		return 0;
	
	return record->from_which_chrom(i, is_mat);
}
