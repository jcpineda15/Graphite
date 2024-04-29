#ifndef __VCFHETEROHOMOPP
#define __VCFHETEROHOMOPP

#include "VCFImpFamily.h"
#include "VCFFillable.h"
#include "group.h"
#include "Map.h"
#include "option.h"
#include "graph.h"
#include "Baum_Welch_with_fixed_Ts.h"

class Map;
class Family;
class KnownFamily;
class PedigreeTable;
class VCFOriginal;


//////////////////// VCFHeteroHomoPP ////////////////////

class VCFHeteroHomoPP : public VCFBase, public VCFSmallBase,
						public VCFFamilyBase, public VCFMeasurable {
public:
	struct ConfigThread {
		const VCFSmall	*orig_vcf;
		const VCFSmall	*merged_vcf;
		const std::vector<const KnownFamily *>&	families;
		const Map& geno_map;
		const std::size_t	first;
		const int	num_threads;
		std::vector<VCFFillable *>&	results;
		
		ConfigThread(const VCFSmall *o, const VCFSmall *m,
					 const std::vector<const KnownFamily *>& fs,
					 const Map& gmap, std::size_t f,
					 int T, std::vector<VCFFillable *>& rs) :
			 						orig_vcf(o), merged_vcf(m),
			 						families(fs), geno_map(gmap),
			 						first(f), num_threads(T), results(rs) { }
		
		std::size_t size() const { return families.size(); }
	};

protected:
	std::vector<VCFFillableRecord *>	records;
	
public:
	VCFHeteroHomoPP(const std::vector<STRVEC>& h, const STRVEC& s,
						std::vector<VCFFillableRecord *> rs, const Map& m);
	~VCFHeteroHomoPP();
	
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
	bool is_mat_hetero() const;
	const std::vector<VCFFillableRecord *>& get_records() const {
		return records;
	}
	
	void impute();
	void clear_records() { records.clear(); }
	void fill();
	
private:
	double record_cM(std::size_t i) const { return cM(records[i]->pos()); }
	std::string make_seq(std::size_t i) const;
	std::string impute_sample_seq(std::size_t i,
								const std::vector<double>& cMs, double min_c);
	std::string update_each(std::size_t i, std::size_t j, char c);
	void update(std::size_t i, const STRVEC& seqs);
	
	void impute_core(const RecordSet& record_set);
	VCFFillableRecord *find_prev_record(
							const Groups *groups, int i, FillType g);
	VCFFillableRecord *find_next_record(
							const Groups *groups, int i, FillType g);
	
public:
	static std::map<FillType, std::vector<VCFFillableRecord *>>
				classify_records(const std::vector<VCFFamilyRecord *>& records);
	static VCFFillable *merge_vcf(const VCFHeteroHomoPP *mat_vcf,
								  const VCFHeteroHomoPP *pat_vcf,
				 const std::vector<VCFFillableRecord *>& homohomo_records,
				 const std::vector<VCFFillableRecord *>& heterohetero_records);
	static VCFFillable *impute_by_parents(const VCFSmallBase *orig_vcf,
										const VCFSmallBase *imputed_vcf,
										const STRVEC& samples, const Map& gmap);
	static std::vector<VCFFillable *> impute_vcfs(
						const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
						const std::vector<const KnownFamily *>& families,
						const Map& geno_map, int num_threads);
	
	static VCFFillableRecord *merge_record(const VCFRecord *record1,
											const VCFRecord *record2,
											const STRVEC& samples, int i,
											const TypeDeterminer *td);
	static VCFFillableRecord *fill_NA(VCFRecord *record1,
											const STRVEC& samples, int i);
	static VCFHeteroHomoPP *merge(const VCFSmallBase *vcf_parents,
									const VCFSmallBase *vcf_progenies,
									const STRVEC& samples,
									const Map& m, const Option *option);
	
	
private:
	static bool is_all_same_without_N(const std::string& seq);
	static std::string create_same_color_string(const std::string& seq);
	static std::pair<ParentComb, FillType> classify_record(
												VCFFamilyRecord *record);
	static void impute_in_thread(void *config);
};
#endif
