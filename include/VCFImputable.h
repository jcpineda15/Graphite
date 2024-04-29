#ifndef __VCFIMPUTABLE
#define __VCFIMPUTABLE

#include "VCF.h"
#include "Map.h"
#include "Haplotype.h"


//////////////////// VCFImputable ////////////////////

class VCFImputable : public VCFSmallBase, public VCFMeasurable {
public:
	VCFImputable(const Map& map_) : VCFSmallBase(), VCFMeasurable(map_) { }
	virtual ~VCFImputable() { };
	
	virtual std::vector<Haplotype> collect_haplotypes_mat(
										std::size_t sample_index) const = 0;
	virtual std::vector<Haplotype> collect_haplotypes_pat(
										std::size_t sample_index) const = 0;
	
	std::vector<int> get_int_gts(std::size_t sample_index) const;
	Haplotype clip_haplotype(std::size_t sample_index, int side) const;
	bool is_block(std::size_t first, std::size_t last) const {
		const double	length = cM(get_record(last)->pos()) -
									cM(get_record(first)->pos());
		if(length < 1.0)
			return true;
		else
			return last - first < 10 and length < 10.0;
	}
	
	template<typename T>
	static std::vector<T *> divide_by_cM(const T *vcf) {
		// If more than 1 cM but not more than 10 pieces and not more than 10 cM
		// are considered to be a block
		std::vector<T *>	vcfs;
		std::size_t	first = 0;
		for(std::size_t i = 1; i < vcf->size(); ++i) {
			if(!vcf->is_block(first, i)) {
				vcfs.push_back(vcf->divide_by_positions(first, i));
				first = i;
			}
		}
		vcfs.push_back(vcf->divide_by_positions(first, vcf->size()));
		return vcfs;
	}
	
	HaplotypePair impute_cM_each_sample(HaplotypePair prev_hap,
										std::size_t sample_index,
										bool exec, bool modify_genotypes);
	// phasing, but not correct
	void set_GT_unmodify(std::size_t i, std::size_t sample_id,
											const std::string& new_GT);
	void set_haplotype(HaplotypePair hap, std::size_t sample_index,
												bool modify_genotypes);
};
#endif
