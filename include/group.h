// group.h
#ifndef __GROUP
#define __GROUP

#include <vector>

#include "VCFImpFamily.h"
#include "RecordSet.h"

class VCFFillableRecord;


//////////////////// Groups ////////////////////

class Groups {
public:
	using Group = std::pair<FillType, std::vector<VCFFillableRecord *>>;
	
private:
	const std::vector<Group>	groups;
	
public:
	Groups(const std::vector<Group>& g) : groups(g) { }
	~Groups() { }
	
	std::size_t size() const { return groups.size(); }
	const Group& get_group(std::size_t i) const { return groups[i]; }
	const FillType get_type(std::size_t i) const { return groups[i].first; }
	const std::vector<VCFFillableRecord *>& get_records(std::size_t i) const {
													return groups[i].second; }
	VCFFillableRecord *find_prev_record(std::size_t i, FillType g) const;
	VCFFillableRecord *find_next_record(std::size_t i, FillType g) const;
	std::vector<RecordSet> create_record_sets() const;
	
	static Groups *create(const std::vector<VCFFillableRecord *>& records);
};
#endif
