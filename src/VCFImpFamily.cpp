#include <map>
#include <set>
#include <cassert>
#include "../include/VCFImpFamily.h"
#include "../include/common.h"

using namespace std;

void VCFImpFamilyRecord::set_00x11_parents(int i, int gt) {
	assert(gt != 1);
	const string	gt1 = gt == 0 ? "0|0" : "1|1";
	const string	gt2 = gt == 0 ? "1|1" : "0|0";
	const int	j = i == 0 ? 1 : 0;
	this->set_GT(i, gt1);
	this->set_GT(j, gt2);
	
	// If the parent's genotypes are swapped, change the progenies' genotypes.
	const string	progeny_GT = (i == 0) ^ (gt == 2) ? "0|1" : "1|0";
	for(int k = 2; k < (int)num_samples(); ++k)
		this->set_GT(k, progeny_GT);
}

int VCFImpFamilyRecord::get_records_type(const vector<RecordWithPos>& records) {
	int	type = -1;
	for(auto q = records.begin(); q != records.end(); ++q) {
		auto	*record = q->first;
		if(record->is_fixed()) {
			if(type == -1)
				type = 0;
			else if(type == 1)
				return 2;
		}
		else {
			if(type == -1)
				type = 1;
			else if(type == 0)
				return 2;
		}
	}
	return type;
}

size_t VCFImpFamilyRecord::find_fixed_index(
				const vector<std::pair<int, vector<RecordWithPos>>>& items) {
	size_t	all_fixed_index = string::npos;
	for(auto p = items.begin(); p != items.end(); ++p) {
		// 0: all fixed 1: all not fixed 2: otherwise
		const int	records_type = get_records_type(p->second);
		if(records_type == 2)
			return string::npos;
		else if(records_type == 0) {
			if(all_fixed_index == string::npos)
				all_fixed_index = p - items.begin();
			else
				return string::npos;
		}
	}
	return all_fixed_index;
}

std::pair<int, vector<VCFImpFamilyRecord::RecordWithPos>>
VCFImpFamilyRecord::which_is_fixed(const vector<RecordWithPos>& v) {
	vector<int>	gts;
	for(auto p = v.begin(); p != v.end(); ++p) {
		gts.push_back(p->first->get_int_gt(p->second));
	}
	if(Common::is_all_same(gts))
		return make_pair(gts.front(), vector<RecordWithPos>());
	
	// collect reads in Genotype
	// ignore ./.
	map<int, vector<RecordWithPos>>	dic;
	for(size_t k = 0U; k < v.size(); ++k) {
		dic[gts[k]].push_back(v[k]);
	}
	const vector<std::pair<int, vector<RecordWithPos>>>
	items(dic.begin(), dic.end());
	
	// What number is all fixed?
	const size_t	fixed_index = find_fixed_index(items);
	if(fixed_index == string::npos)
		return make_pair(gts.front(), vector<RecordWithPos>());
	
	const int	fixed_GT = items[fixed_index].first;
	if(fixed_GT == 1)	// 0/1
		return make_pair(gts.front(), vector<RecordWithPos>());
	
	vector<RecordWithPos>	wrongs;
	for(size_t i = 0U; i < items.size(); ++i) {
		if(i == (size_t)fixed_index)
			continue;
		const auto&	records = items[i].second;
		for(auto p = records.begin(); p != records.end(); ++p) {
			if(p->first->is_00x11())
				wrongs.push_back(*p);
		}
	}
	return make_pair(fixed_GT, wrongs);
}

void VCFImpFamilyRecord::modify_00x11(
						const vector<VCFImpFamilyRecord *>& records) {
	// collect records by sample
	map<string, vector<std::pair<VCFImpFamilyRecord *, int>>>	dic;
	for(auto p = records.begin(); p != records.end(); ++p) {
		dic[(*p)->get_samples()[0]].push_back(make_pair(*p, 0));
		dic[(*p)->get_samples()[1]].push_back(make_pair(*p, 1));
	}
	
	for(auto p = dic.begin(); p != dic.end(); ++p) {
		if(p->second.size() < 2)
			continue;
		
		auto	v = which_is_fixed(p->second);
		const int	fixed_gt = v.first;
		const auto	wrongs = v.second;
		for(auto q = wrongs.begin(); q != wrongs.end(); ++q) {
			auto	*record = q->first;
			int		i = q->second;
			record->set_00x11_parents(i, fixed_gt);
		}
	}
}
