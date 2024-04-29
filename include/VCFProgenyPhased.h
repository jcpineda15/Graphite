#ifndef __VCFPROGENYPHASED
#define __VCFPROGENYPHASED

#include "VCFFamily.h"
#include "VCFImputable.h"

class Family;
class KnownFamily;


//////////////////// VCFOneParentPhased ////////////////////

class VCFProgenyPhased : public VCFBase,
							public VCFFamilyBase, public VCFImputable {
public:
	typedef std::pair<const KnownFamily *, std::size_t>	PAIR;
	
	struct ConfigThread {
		const VCFSmall	*orig_vcf;
		const VCFSmall	*merged_vcf;
		const std::vector<PAIR>&	families;
		const Map& map;
		const VCFSmall *ref_vcf;
		const std::size_t	first;
		const int	num_threads;
		std::vector<VCFProgenyPhased *>&	results;
		
		ConfigThread(const VCFSmall *o, const VCFSmall *m,
					 const std::vector<PAIR>& fs,
					 const Map& map_, const VCFSmall *ref,
					 std::size_t f, int T,
					 std::vector<VCFProgenyPhased *>& rs) :
			 						orig_vcf(o), merged_vcf(m),
			 						families(fs), map(map_), ref_vcf(ref),
			 						first(f), num_threads(T), results(rs) { }
		
		std::size_t size() const { return families.size(); }
	};
	
protected:
	std::vector<VCFFamilyRecord *>	records;
	int	selection;
	const VCFSmall	*ref_vcf;
	
public:
	VCFProgenyPhased(const std::vector<STRVEC>& h, const STRVEC& s,
									std::vector<VCFFamilyRecord *> rs,
									const Map& gmap, const VCFSmall *ref);
	~VCFProgenyPhased();
	
	///// virtual methods for VFSmallBase /////
	const std::vector<STRVEC>& get_header() const {
		return VCFBase::get_header();
	}
	const STRVEC& get_samples() const { return VCFBase::get_samples(); }
	std::size_t size() const { return records.size(); }
	VCFRecord *get_record(std::size_t i) const {
		return records[i];
	}
	
	///// virtual methods for VFSmallBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const {
		return records[i];
	}
	
	///// virtual methods for VCFImputable /////
	std::vector<Haplotype> collect_haplotypes_mat(
									std::size_t sample_index) const;
	std::vector<Haplotype> collect_haplotypes_pat(
									std::size_t sample_index) const;
	void set_gts(const std::vector<std::string>& gts, std::size_t sample_index);
	
	///// non-virtual methods /////
	VCFProgenyPhased *divide_by_positions(std::size_t first,
											std::size_t last) const;
	void impute();
	
private:
	Haplotype clip_ref_haplotype(std::size_t sample_index, int side) const;
	std::vector<Haplotype> collect_haplotypes_from_phased_progeny(
															int side) const;
	std::vector<Haplotype> collect_haplotype_from_refs() const;
	std::size_t collect_known_parents_indices() const;
	std::vector<HaplotypePair> impute_core(
								const std::vector<VCFProgenyPhased *>& vcf_cMs);
	int score_whole(const std::vector<HaplotypePair>& haps,
						std::size_t sample_index,
						const std::vector<VCFProgenyPhased *>& vcf_cMs) const;
	int sum_score(const std::vector<HaplotypePair>& hap_pairs,
						const std::vector<VCFProgenyPhased *>& vcf_cMs,
						const std::size_t sample_index) const;
	
	void inverse_selection() {
		selection = selection == 0 ? 1 : 0;
	}
	
private:
	static void impute_in_thread(void *config);
	
public:
	static VCFProgenyPhased *impute_by_progeny(const VCFSmall *orig_vcf,
												const VCFSmall *imputed_vcf,
												const STRVEC& samples,
												const std::size_t& ppi,
												const Map& map_,
												const VCFSmall *ref_vcf);
	static std::vector<VCFProgenyPhased *> impute_all_by_progeny(
							const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
							const std::vector<PAIR>& families, const Map& gmap,
							const VCFSmall *ref_vcf, int num_threads);
	static VCFProgenyPhased *create(const VCFSmall *orig_vcf,
							const VCFSmall *merged_vcf,
							const STRVEC& family_samples, size_t phased_index,
							const Map& map_, const VCFSmall *ref_vcf);
};
#endif
