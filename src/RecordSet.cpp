#include <iostream>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "../include/RecordSet.h"
#include "../include/common.h"

using namespace std;


//////////////////// RecordSet ////////////////////

int RecordSet::prev_mat_from(std::size_t i) const {
	if(prev_mat_record == NULL)
		return 0;
	return prev_mat_record->from_which_chrom(i, true);
}
int RecordSet::next_mat_from(std::size_t i) const {
	if(next_mat_record == NULL)
		return 0;
	return next_mat_record->from_which_chrom(i, true);
}
int RecordSet::prev_pat_from(std::size_t i) const {
	if(prev_pat_record == NULL)
		return 0;
	return prev_pat_record->from_which_chrom(i, false);
}
int RecordSet::next_pat_from(std::size_t i) const {
	if(next_pat_record == NULL)
		return 0;
	return next_pat_record->from_which_chrom(i, false);
}

bool RecordSet::is_mat_prev_near() const {
	return record->pos() * 2 < prev_mat_record->pos() + next_mat_record->pos();
}

bool RecordSet::is_pat_prev_near() const {
	return record->pos() * 2 < prev_pat_record->pos() + next_pat_record->pos();
}

int RecordSet::near_mat_from(size_t i) const {
	return is_mat_prev_near() ? prev_mat_from(i) : next_mat_from(i);
}

int RecordSet::near_pat_from(size_t i) const {
	return is_pat_prev_near() ? prev_pat_from(i) : next_pat_from(i);
}

RecordSet::Pair RecordSet::select_nearest_froms(
								const vector<Pair>& pairs, size_t i) const {
	if(pairs.size() == 4U) {
		return Pair(near_mat_from(i), near_pat_from(i));
	}
	else if(pairs[0].first == pairs[1].first) {			// same mat
		if(is_pat_prev_near())
			return Pair(pairs[0].first, prev_pat_from(i));
		else
			return Pair(pairs[0].first, next_pat_from(i));
	}
	else if(pairs[0].second == pairs[1].second) {	// same pat
		if(is_mat_prev_near())
			return Pair(prev_mat_from(i), pairs[0].second);
		else
			return Pair(next_mat_from(i), pairs[1].second);
	}
	else {	// Both parents crossover (rarely)
		return Pair(near_mat_from(i), near_pat_from(i));
	}
}


// Select a pair from both parent's Haplotype
// Prefer pairs that do not change the genotype
RecordSet::Pair RecordSet::select_pair(const vector<Pair>& pairs,
												size_t i, bool selected) const {
	if(pairs.empty())
		return Pair(0, 0);
	else if(pairs.size() == 1U)
		return pairs.front();
	else if(!Genotype::is_valid(record->get_gt(i), record->mat_int_gt(),
													record->pat_int_gt()))
		return select_nearest_froms(pairs, i);
	else if(selected)
		return select_nearest_froms(pairs, i);
	
	// Collect pairs that are identical to the original Genotype
	vector<Pair>	new_pairs;
	for(auto p = pairs.begin(); p != pairs.end(); ++p) {
		string	parent_gt = record->gt_from_parent(p->first, p->second);
		if(Genotype::sum_gt(record->get_gt(i)) == Genotype::sum_gt(parent_gt))
			new_pairs.push_back(*p);
	}
	const Pair	selected_pair = select_pair(new_pairs, i, true);
	if(selected_pair.first != 0)
		return selected_pair;
	else
		return select_pair(pairs, i, true);
}

// Likelihood based on which of the parent Haplotypes
// came from before and after the record
vector<double> RecordSet::likelihoods_from_which_chrom(
										int prev_from, int next_from) const {
	static const double ps[] = {
		0.5, 0.9, 0.1, 0.9, 0.99, 0.5, 0.1, 0.5, 0.01
	};
	const double	prob1 = ps[prev_from+next_from*3];
	vector<double>	probs = { prob1, 1.0 - prob1 };
	return probs;
}

vector<double> RecordSet::likelihoods_from_which_chrom(
												size_t i, bool is_mat) const {
	if(is_mat)
		return likelihoods_from_which_chrom(prev_mat_from(i), next_mat_from(i));
	else
		return likelihoods_from_which_chrom(prev_pat_from(i), next_pat_from(i));
}

double RecordSet::likelihood_each(const string& gt,
									const vector<double>& probs_mat,
									const vector<double>& probs_pat,
									int mat_phasing, int pat_phasing) const {
	const int	sum = Genotype::sum_gt(gt);
	double	likelihood = 0.0;
	for(int k = 0; k < 4; ++k) {
		const int	i = k >> 1;
		const int	j = k & 1;
		// does this pair match genotype?
		if((((mat_phasing >> i) & 1) + ((pat_phasing >> j) & 1)) == sum)
			likelihood += probs_mat[i] * probs_pat[j];
	}
	if(likelihood == 0.0)	// no matching pair
		return log(0.0001);
	else
		return log(likelihood);
}

double RecordSet::compute_phasing_likelihood_each(size_t i,
									int mat_phasing, int pat_phasing) const {
	const auto	probs_mat = this->likelihoods_from_which_chrom(i, true);
	const auto	probs_pat = this->likelihoods_from_which_chrom(i, false);
	return this->likelihood_each(this->gt(i), probs_mat, probs_pat,
											mat_phasing, pat_phasing);
}

double RecordSet::compute_phasing_likelihood(int mat_phasing,
														int pat_phasing) const {
	double	ll = 0.0;
	for(int i = 2; i < (int)record->num_samples(); ++i) {
		if(record->get_GT(i) == "./.")
			ll += log(0.0001);	// to match the Python version
		else
			ll += compute_phasing_likelihood_each(i, mat_phasing, pat_phasing);
	}
	return ll;
}

bool RecordSet::is_prev_nearer(bool is_mat) const {
	if(record == NULL)
		return false;
	else if(is_mat) {
		// どちらかがNULLの場合、ここには来ないので、適当に処理する
		if(prev_mat_record == NULL || next_mat_record == NULL)
			return false;
		else
			return is_mat_prev_near();
	}
	else {
		if(prev_pat_record == NULL || next_pat_record == NULL)
			return false;
		else
			return is_pat_prev_near();
	}
}

pair<int, int> RecordSet::select_phasing(
							const vector<pair<int, int>>& candidates) const {
	if(candidates.size() == 1)
		return candidates.front();
	
	const int	mat_int_gt = record->mat_int_gt();
	const int	pat_int_gt = record->pat_int_gt();
	
	vector<tuple<int, int, int>>	v;	// [(score, mat_phasing, pat_phasing)]
	for(auto p = candidates.begin(); p != candidates.end(); ++p) {
		const int	mat_int_gt1 = (p->first >> 1) + (p->first & 1);
		const int	pat_int_gt1 = (p->second >> 1) + (p->second & 1);
		const int	score = abs(mat_int_gt1 - mat_int_gt) +
							abs(pat_int_gt1 - pat_int_gt);
		v.push_back(make_tuple(score, p->first, p->second));
	}
	// If the score is the same, phasing is chosen from the smaller.
	// want to choice at random
	const auto	p = std::min_element(v.begin(), v.end());
	return make_pair(get<1>(*p), get<2>(*p));
}

pair<int, int> RecordSet::determine_phasing_core(
							const vector<tuple<double, int, int>>& lls) const {
	// collect max or near likelihoods
	// for numerical error
	vector<pair<int, int>>	candidates;
	const double	max_ll = get<0>(lls.back());
	for(size_t i = lls.size() - 1; i < lls.size(); --i) {
		const double	ll = get<0>(lls[i]);
		if(ll > max_ll - 1e-9)
			candidates.push_back(make_pair(get<1>(lls[i]), get<2>(lls[i])));
		else
			break;
	}
	
	return this->select_phasing(candidates);
}

void RecordSet::determine_phasing() const {
	if(this->record == NULL)
		return;
	
	const vector<pair<int, int>>	phasing = record->possible_phasings();
	vector<tuple<double, int, int>>	lls;
	for(auto p = phasing.begin(); p != phasing.end(); ++p) {
		const double	ll = compute_phasing_likelihood(p->first, p->second);
		lls.push_back(make_tuple(ll, p->first, p->second));
	}
	
	// Measures for when there is almost the same ll as the largest ll
	std::sort(lls.begin(), lls.end());
	const auto	p = this->determine_phasing_core(lls);
	const int	mat_phasing = p.first;
	const int	pat_phasing = p.second;
	
	// is it OK?
	static const string	gts[] = { "0|0", "1|0", "0|1", "1|1" };
	record->set_mat_GT(gts[mat_phasing]);
	record->set_pat_GT(gts[pat_phasing]);
}

int RecordSet::select_from(int from1, int from2,
										const VCFRecord *record1,
										const VCFRecord *record2) const {
	// assume either is not zero
	// That is, there is a record on either side.
	if(from1 == 0)
		return from2;
	else if(from2 == 0)
		return from1;
	else {	// there are records on both sides
		// select a close record
		if(record->pos() * 2 < record1->pos() + record2->pos())
			return from1;
		else
			return from2;
	}
}

string RecordSet::modify_gt(size_t i) const {
	const int	prev_mat_from = this->prev_mat_from(i);
	const int	next_mat_from = this->next_mat_from(i);
	const int	prev_pat_from = this->prev_pat_from(i);
	const int	next_pat_from = this->next_pat_from(i);
	if((prev_mat_from == 0 && next_mat_from == 0) ||
							(prev_pat_from == 0 && next_pat_from == 0)) {
		return record->get_gt(i);
	}
	
	vector<pair<int,int>>	pairs_ = {
				pair<int,int>(prev_mat_from, prev_pat_from),
				pair<int,int>(prev_mat_from, next_pat_from),
				pair<int,int>(next_mat_from, prev_pat_from),
				pair<int,int>(next_mat_from, next_pat_from)
	};
	vector<pair<int,int>>	pairs__;
	for(auto p = pairs_.begin(); p != pairs_.end(); ++p) {
		if(p->first != 0 && p->second != 0)
			pairs__.push_back(*p);
	}
	const auto	pairs = Common::unique_vector(pairs__);
	const auto	p = select_pair(pairs, i);
	const int	mat_from = p.first;
	const int	pat_from = p.second;
	if(mat_from == 0 && pat_from == 0) {
		if(prev_mat_from != 0 || next_mat_from != 0) {
			const int	mat_from_selected = select_from(
											prev_mat_from, next_mat_from,
											prev_mat_record, next_mat_record);
			return record->gt_from_mat(mat_from_selected, i + 9);
		}
		else if(prev_pat_from != 0 || next_pat_from != 0) {
			const int	pat_from_selected = select_from(
											prev_pat_from, next_pat_from,
											prev_pat_record, next_pat_record);
			return record->gt_from_pat(pat_from_selected, i + 9);
		}
		else {
			// if both are 0, do nothing
			return record->get_GT(i);
		}
	}
	else {
		return record->gt_from_parent(mat_from, pat_from);
	}
}

void RecordSet::impute(bool necessary_parents_phasing) const {
	if(necessary_parents_phasing)
		this->determine_phasing();
	this->impute_core();
}

void RecordSet::impute_core() const {
	STRVEC	new_gts;
	for(size_t i = 2; i < record->num_samples(); ++i) {
		new_gts.push_back(modify_gt(i));
	}
	record->modify_gts(new_gts);
	record->modify_parents_type();
}

int RecordSet::select_mat(const vector<Pair>& pairs) const {
	if(pairs.size() == 1)
		return pairs.front().first;
	else if(this->prev_mat_record == NULL || this->next_mat_record == NULL)
		return 0;
	else if(this->is_mat_prev_near())
		return pairs.front().first;		// prev_mat_from
	else
		return pairs.back().first;		// next_mat_from
}

void RecordSet::impute_NA_mat_each(size_t i) const {
	const int	prev_mat_from = VCFFillableRecord::from_which_chrom_mat(
													this->prev_mat_record, i);
	const int	next_mat_from = VCFFillableRecord::from_which_chrom_mat(
													this->next_mat_record, i);
	
	// pat can be 1 or 2 because pat is homozygous
	vector<Pair>	pairs;	// [(mat_from, pat_from)]
	if(this->prev_mat_record != NULL)
		pairs.push_back(Pair(prev_mat_from, 1));
	if(this->next_mat_record != NULL && next_mat_from != prev_mat_from)
		pairs.push_back(Pair(next_mat_from, 1));
	if(pairs.empty())
		return;
	
	const int	mat_from = this->select_mat(pairs);
	this->record->set_GT(i, this->record->gt_from_parent(mat_from, 1));
}

int RecordSet::select_pat(const vector<Pair>& pairs) const {
	if(pairs.size() == 1)
		return pairs.front().second;
	else if(this->prev_pat_record == NULL || this->next_pat_record == NULL)
		return 0;
	else if(this->is_pat_prev_near())
		return pairs.front().second;	// prev_pat_from
	else
		return pairs.back().second;		// next_pat_from
}

void RecordSet::impute_NA_pat_each(size_t i) const {
	const int	prev_pat_from = VCFFillableRecord::from_which_chrom_pat(
													this->prev_pat_record, i);
	const int	next_pat_from = VCFFillableRecord::from_which_chrom_pat(
													this->next_pat_record, i);
	
	vector<Pair>	pairs;
	if(this->prev_pat_record != NULL)
		pairs.push_back(Pair(1, prev_pat_from));
	if(this->next_pat_record != NULL && next_pat_from != prev_pat_from)
		pairs.push_back(Pair(1, next_pat_from));
	if(pairs.empty())
		return;
	
	const int	pat_from = this->select_pat(pairs);
	this->record->set_GT(i, this->record->gt_from_parent(1, pat_from));
}

int RecordSet::determine_mat_from(size_t i) const {
	const string	mat_gt1 = gt_each(i, prev_mat_record);
	const string	mat_gt2 = gt_each(i, next_mat_record);
	const int	prev_mat_from = from_which_chrom_prev_mat(mat_gt1);
	const int	next_mat_from = from_which_chrom_next_mat(mat_gt2);
	// とりあえず、両側Noneはないと仮定する
	if(prev_mat_from == 0)
		return next_mat_from;
	else if(next_mat_from == 0)
		return prev_mat_from;
	else if(prev_mat_from == next_mat_from)
		return prev_mat_from;
	else if(is_prev_nearer(true))
		return prev_mat_from;
	else
		return next_mat_from;
}

int RecordSet::determine_pat_from(size_t i) const {
	const string	pat_gt1 = gt_each(i, prev_pat_record);
	const string	pat_gt2 = gt_each(i, next_pat_record);
	const int	prev_pat_from = from_which_chrom_prev_pat(pat_gt1);
	const int	next_pat_from = from_which_chrom_next_pat(pat_gt2);
	// とりあえず、両側Noneはないと仮定する
	if(prev_pat_from == 0)
		return next_pat_from;
	else if(next_pat_from == 0)
		return prev_pat_from;
	else if(prev_pat_from == next_pat_from)
		return prev_pat_from;
	else if(is_prev_nearer(false))
		return prev_pat_from;
	else
		return next_pat_from;
}
