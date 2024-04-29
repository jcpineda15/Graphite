// group.cpp
#include "../include/group.h"
#include "../include/VCFFillableRecord.h"

using namespace std;


//////////////////// Groups ////////////////////

VCFFillableRecord *Groups::find_prev_record(size_t i, FillType g) const {
	if(i == 0)
		return NULL;
	
	for(size_t j = i-1; j < i; --j) {	// underflow
		const FillType	key = get_group(j).first;
		const auto&	records = get_group(j).second;
		if(key == g)
			return records.back();
	}
	return NULL;
}

VCFFillableRecord *Groups::find_next_record(size_t i, FillType g) const {
	for(size_t j = i+1; j < size(); ++j) {
		const FillType	key = get_group(j).first;
		const auto&	records = get_group(j).second;
		if(key == g)
			return records.front();
	}
	return NULL;
}

vector<RecordSet> Groups::create_record_sets() const {
	vector<RecordSet>	record_sets;
	for(size_t i = 0; i < size(); ++i) {
		const FillType	key = get_group(i).first;
		const auto&	records = get_group(i).second;
		if(key == FillType::MAT || key == FillType::PAT)
			continue;
		auto	*prev_mat_record = find_prev_record(i, FillType::MAT);
		auto	*next_mat_record = find_next_record(i, FillType::MAT);
		auto	*prev_pat_record = find_prev_record(i, FillType::PAT);
		auto	*next_pat_record = find_next_record(i, FillType::PAT);
		for(auto p = records.begin(); p != records.end(); ++p) {
			record_sets.push_back(RecordSet(*p, prev_mat_record,
						next_mat_record, prev_pat_record, next_pat_record));
		}
	}
	return record_sets;
}

Groups *Groups::create(const vector<VCFFillableRecord *>& records) {
	vector<Group>	groups;
	auto	current_type = records.front()->get_type();
	vector<VCFFillableRecord *>	group(1U, records.front());
	for(auto p = records.begin() + 1; p != records.end(); ++p) {
		auto	*record = *p;
		if(record->get_type() != current_type) {
			groups.push_back(Group(current_type, group));
			group.clear();
			current_type = record->get_type();
		}
		group.push_back(record);
	}
	groups.push_back(Group(current_type, group));
	return new Groups(groups);
}
