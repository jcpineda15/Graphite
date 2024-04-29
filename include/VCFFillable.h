#ifndef __VCFFILLABLE
#define __VCFFILLABLE

#include "RecordSet.h"

class VCFHeteroHomo;
class Option;


//////////////////// VCFFillable ////////////////////

class VCFFillable : public VCFBase, public VCFSmallBase, public VCFFamilyBase {
	using PosWithChr = std::tuple<int,ll,std::string>;
	
	
public:
	using Position = std::tuple<int, ll, std::string>;
	using Item = std::pair<std::vector<VCFHeteroHomo *>,
							std::vector<VCFImpFamilyRecord *>>;
	using ImpRecords = std::map<Parents, std::vector<VCFImpFamilyRecord *>>;
	
	struct ConfigThreadPhase {
		const std::size_t	first;
		const std::size_t	num_threads;
		const std::vector<RecordSet>&	record_sets;
		VCFFillable	*vcf;
		
		ConfigThreadPhase(std::size_t f, std::size_t n,
							const std::vector<RecordSet>& r, VCFFillable *v) :
							first(f), num_threads(n), record_sets(r), vcf(v) { }
	};
	
public:
	
	std::vector<VCFFillableRecord *>	records;
	
public:
	VCFFillable(const std::vector<STRVEC>& h, const STRVEC& s,
									std::vector<VCFFillableRecord *> rs);
	~VCFFillable();
	
	///// virtual methods for VCFSmallBase /////
	std::size_t size() const { return records.size(); }
	VCFRecord *get_record(std::size_t i) const {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	const std::vector<STRVEC>& get_header() const {
		return VCFBase::get_header();
	}
	const STRVEC& get_samples() const { return VCFBase::get_samples(); }
	VCFFamilyRecord *get_family_record(std::size_t i) const {
		return records[i];
	}
	
	///// non-virtual methods /////
	VCFFillableRecord *get_fillable_record(std::size_t i) const {
		return records[i];
	}
	const std::vector<VCFFillableRecord *>& get_records() const {
		return records;
	}
	
	void modify(int num_threads);
	void phase_hetero_hetero();
	VCFFillable *create_from_header() const;
	void set_records(const std::vector<VCFFillableRecord *>& rs) {
		records = rs;
	}
	void clear_records() { records.clear(); }
	
private:
	template<typename Iter>
	VCFFillableRecord *find_neighbor_same_type_record(
					std::size_t i, std::size_t c, Iter first, Iter last) const;
	VCFFillableRecord *find_prev_same_type_record(std::size_t i,
													std::size_t c) const;
	VCFFillableRecord *find_next_same_type_record(std::size_t i,
													std::size_t c) const;
	
	const RecordSet *create_recordset(std::size_t i,
										std::size_t c, bool is_mat) const;
	void impute_NA_mat_each(std::size_t i, std::size_t c);
	void impute_NA_mat(std::size_t i);
	
	void impute_NA_pat_each(std::size_t i, std::size_t c);
	void impute_NA_pat(std::size_t i);
	
	std::pair<int, int> find_prev_mat_from(int i, int c) const;
	std::pair<int, int> find_next_mat_from(int i, int c) const;
	std::pair<int, int> find_prev_pat_from(int i, int c) const;
	std::pair<int, int> find_next_pat_from(int i, int c) const;
	int select_from(const std::pair<int, int>& f1,
					const std::pair<int, int>& f2, int i) const;
	int find_mat_from(int i, int c) const;
	int find_pat_from(int i, int c) const;
	void impute_others(int i);
	
public:
	static VCFSmall *merge(const std::vector<VCFFillable *>& vcfs,
											const STRVEC& orig_samples);
	static std::vector<VCFFillable *> fill_all(
				std::map<Parents, std::vector<VCFHeteroHomo *>>& imputed_vcfs,
				ImpRecords& other_records, int num_threads);
	static VCFFillable *fill(const std::vector<VCFHeteroHomo *>& vcfs,
							 const std::vector<VCFImpFamilyRecord *>& records,
							 int num_threads);
	static std::vector<VCFFillableRecord *> merge_records(
							const std::vector<VCFHeteroHomo *>& vcfs,
							const std::vector<VCFImpFamilyRecord *>& records,
							bool all_out);
	
private:
	static std::vector<std::vector<VCFFillableRecord *>>
		collect_records(const std::vector<VCFFillable *>& vcfs);
	static std::pair<STRVEC, std::vector<std::vector<std::pair<int, int>>>>
			integrate_samples(const std::vector<STRVEC>& sss,
											const STRVEC& orig_samples);
	// Integrate the VCF so that duplicate samples are one
	static VCFSmall *integrate(const VCFFillable *vcf,
					const std::vector<std::vector<VCFFillableRecord *>>& rss,
					const STRVEC& orig_samples);
	static std::vector<VCFFillable *> fill_parellel(std::vector<Item>& items,
															int num_threads);
	static void delete_items(const std::vector<Item>& items);
	static void phase_in_thread(void *config);
	static void fill_in_thread(void *config);
};


//////////////////// ConfigFillThread ////////////////////

struct ConfigFillThread {
	const std::vector<VCFFillable::Item>&	items;
	const bool	all_out;
	const std::size_t	first;
	const int	num_threads;
	std::vector<VCFFillable *>&	filled_vcfs;
	
	ConfigFillThread(const std::vector<VCFFillable::Item>& items_,
						bool ao, int f, int n,
						std::vector<VCFFillable *>& results) :
									items(items_), all_out(ao), first(f),
									num_threads(n), filled_vcfs(results) { }
	
	std::size_t size() const { return items.size(); }
};
#endif
