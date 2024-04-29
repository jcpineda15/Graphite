#ifndef __IMPUTEPROGONLY
#define __IMPUTEPROGONLY

#include <vector>
#include <string>
#include <utility>

#include "VCFImpFamily.h"

class TypeDeterminer;
class VCFHeteroHomoRecord;
class VCFImpFamilyRecord;
class VCFFamily;
class VCFHeteroHomoPP;
class VCFFillableRecord;
class Map;
class Option;


//////////////////// ImputeProgOnly ////////////////////

namespace ImputeProgOnly {
	void classify_record(std::size_t i, VCFFamilyRecord *record,
							const STRVEC& samples, const TypeDeterminer *td,
							std::vector<VCFHeteroHomoRecord *>& heho_records,
							std::vector<VCFImpFamilyRecord *>& other_records);
	std::pair<std::vector<VCFHeteroHomoRecord *>,
			  std::vector<VCFImpFamilyRecord *>>
	classify_records(VCFFamily *vcf, const STRVEC& samples, Option *option);
	VCFHeteroHomoPP *merge_vcf(
					std::map<FillType, std::vector<VCFFillableRecord *>>& rss,
					const std::vector<STRVEC>& header,
					const STRVEC& samples, const Map& gmap);
	VCFHeteroHomoPP *impute_prog_vcf_chr(const VCFSmallBase *parent_vcf,
										const VCFSmallBase *prog_vcf,
										const Map& gmap, const Option *Option);
}
#endif
