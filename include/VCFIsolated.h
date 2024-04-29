#ifndef __VCFISOLATED
#define __VCFISOLATED

#include "VCF.h"
#include "Map.h"
#include "VCFImputable.h"


//////////////////// VCFIsolated ////////////////////

// imputed samples at the beginning, followed by reference samples
class VCFIsolated : public VCFBase, public VCFImputable {
public:
	struct ConfigThread {
		const std::vector<VCFIsolated *>&	vcfs;
		const std::size_t	first;
		const int	num_threads;
		
		ConfigThread(const std::vector<VCFIsolated *>& vcfs_,
											std::size_t f, int n) :
									vcfs(vcfs_), first(f), num_threads(n) { }
		
		std::size_t size() const { return vcfs.size(); }
	};
	
private:
	std::vector<VCFRecord *>	records;
	const std::size_t	num_imputed_samples;
	const bool	modify_genotypes;
	
public:
	VCFIsolated(const std::vector<STRVEC>& h, const STRVEC& s,
				std::vector<VCFRecord *> rs, std::size_t nis,
				const Map& m, bool modify_genotypes);
	~VCFIsolated();
	
	///// virtual methods for VFSmallBase /////
	const std::vector<STRVEC>& get_header() const {
		return VCFBase::get_header();
	}
	const STRVEC& get_samples() const { return VCFBase::get_samples(); }
	std::size_t size() const { return records.size(); }
	VCFRecord *get_record(std::size_t i) const {
		return records[i];
	}
	
	///// virtual methods for VCFImputable /////
	std::vector<Haplotype> collect_haplotypes_mat(
									std::size_t sample_index) const;
	std::vector<Haplotype> collect_haplotypes_pat(
									std::size_t sample_index) const;
	void set_gts(const std::vector<std::string>& gts, std::size_t sample_index);
	
	///// non-virtual methods /////
	VCFIsolated *divide_by_positions(std::size_t first, std::size_t last) const;
	void impute();
	VCFSmall *extract_isolated_samples() const;
	void add_record(VCFRecord *record) { records.push_back(record); }
	
	std::vector<Haplotype> collect_haplotype_from_refs() const;
	
private:
	std::vector<HaplotypePair> impute_cM(
								const std::vector<HaplotypePair>& prev_haps);
	
public:
	static std::vector<VCFIsolated *> create(const VCFSmall *orig_vcf,
					const VCFSmall *merged_vcf,
					const STRVEC& samples, const STRVEC& references,
					const Map& gmap, bool modify_genotypes, int num_threads);
	static std::vector<VCFSmallBase *> impute_all(
										const std::vector<VCFIsolated *>& vcfs,
										int num_threads);
	
private:
	static std::vector<std::vector<std::size_t>>
			divide_columns(const std::vector<std::size_t>& cs, int num);
	static void impute_parellel(const std::vector<VCFIsolated *>& vcfs,
															int num_threads);
};
#endif
