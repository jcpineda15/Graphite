#ifndef __LARGEFAMILY
#define __LARGEFAMILY

#include <vector>
#include "VCF.h"
#include "ClassifyRecord.h"

class VCFSmall;
class VCFFamily;
class VCFHeteroHomo;
class VCFFillable;
class VCFHeteroHomoRecord;
class VCFImpFamilyRecord;
class Family;
class KnownFamily;
class Map;
class Option;
class TypeDeterminer;

namespace LargeFamily {
	typedef std::pair<std::string,std::string>	Parents;
	typedef std::map<Parents, std::vector<VCFHeteroHomoRecord *>>	HeHoRecords;
	typedef std::map<Parents, std::vector<VCFImpFamilyRecord *>>	ImpRecords;
	
	struct ConfigThreadClassify {
		std::size_t	first;
		std::size_t	num_threads;
		const VCFFamily	*vcf;
		const TypeDeterminer	*td;
		const KnownFamily	*family;
		std::vector<VCFHeteroHomoRecord *>&	heho_records;
		std::vector<VCFImpFamilyRecord *>&	other_records;
		
		ConfigThreadClassify(std::size_t i, std::size_t	n,
							 const VCFFamily *vcf_, const TypeDeterminer *td_,
							 const KnownFamily *fam,
							 std::vector<VCFHeteroHomoRecord *>& heho_rs,
							 std::vector<VCFImpFamilyRecord *>& other_rs) :
					first(i), num_threads(n), vcf(vcf_), td(td_), family(fam),
					heho_records(heho_rs), other_records(other_rs) { }
	};
	
	VCFHeteroHomoRecord *create_heterohomo_record(const STRVEC& v,
									const KnownFamily *family, std::size_t i,
									WrongType wrong_type, ParentComb pc);
	// separate hetero x homo records from others
	void classify_record(std::size_t i, const VCFFamily *vcf,
						 const TypeDeterminer *td, const KnownFamily *family,
						 std::vector<VCFHeteroHomoRecord *>& heho_records,
						 std::vector<VCFImpFamilyRecord *>& other_records);
	
	void classify_records_in_thread(void *config);
	
	// Collect only hetero × homo to impute.
	// Homo × Homo is imputed for each record.
	// The others are imputed by VCFFillable.
	std::pair<std::vector<VCFHeteroHomoRecord *>,
			  std::vector<VCFImpFamilyRecord *>>
	classify_records(const VCFFamily *vcf,
						const KnownFamily *family, const Option *option);
	
	std::vector<VCFFamily *> create_family_vcfs(
							const VCFSmall *orig_vcf,
							const std::vector<const KnownFamily *>& families,
							int num_threads);
	
	void create_vcf_in_thread(void *config);
	
	struct ConfigThreadCreate {
		std::size_t	first;
		std::size_t	num_threads;
		const VCFSmall	*orig_vcf;
		const std::vector<const KnownFamily *>&	families;
		std::vector<VCFFamily *>&	results;
		
		ConfigThreadCreate(int i, int n, const VCFSmall *orig,
							const std::vector<const KnownFamily *>& fams,
							std::vector<VCFFamily *>& res) :
								first(i), num_threads(n), orig_vcf(orig),
								families(fams), results(res) { }
	};
	
	void clean_in_thread(void *config);
	
	struct ConfigThreadClean {
		std::size_t	first;
		std::size_t	num_threads;
		const std::vector<std::vector<VCFHeteroHomoRecord *>>&	recordss;
		const VCFSmall	*orig_vcf;
		const std::vector<const KnownFamily *>&	families;
		const Map&	geno_map;
		std::vector<std::vector<VCFHeteroHomo *>>&	vcfss;
		std::vector<std::vector<VCFHeteroHomoRecord *>>&	unused_recordss;
		
		ConfigThreadClean(int i, int n,
					const std::vector<std::vector<VCFHeteroHomoRecord *>>& rss,
					const VCFSmall *orig,
					const std::vector<const KnownFamily *>& fams, const Map& m,
					std::vector<std::vector<VCFHeteroHomo *>>& vcfss_,
					std::vector<std::vector<VCFHeteroHomoRecord *>>& un) :
								first(i), num_threads(n), recordss(rss),
								orig_vcf(orig), families(fams), geno_map(m),
								vcfss(vcfss_), unused_recordss(un) { }
	};
	
	void fill_in_thread(void *config);
	
	struct ConfigThreadFill {
		std::size_t	first;
		std::size_t	num_threads;
		const std::vector<std::vector<VCFHeteroHomo *>>&	vcfss_heho;
		const std::vector<std::vector<VCFImpFamilyRecord *>>&	other_recordss;
		const int	num_threads_in_family;
		std::vector<VCFFillable *>&	results;
		
		ConfigThreadFill(int i, int n,
				const std::vector<std::vector<VCFHeteroHomo *>>& heho,
				const std::vector<std::vector<VCFImpFamilyRecord *>>&	other,
				int n_f, std::vector<VCFFillable *>& res) :
								first(i), num_threads(n), vcfss_heho(heho),
								other_recordss(other),
								num_threads_in_family(n_f), results(res) { }
	};
	
	std::vector<std::pair<std::string,
						  std::vector<std::pair<std::size_t, std::size_t>>>>
	collect_same_parent_families(
							const std::vector<const KnownFamily *>& families);
	std::vector<std::vector<VCFImpFamilyRecord *>>
	collect_same_position_records(
				const std::vector<std::tuple<std::vector<VCFHeteroHomoRecord *>,
											 std::vector<VCFImpFamilyRecord *>,
											 int>>& rs);
	void modify_00x11_each(
				const std::vector<std::tuple<std::vector<VCFHeteroHomoRecord *>,
											 std::vector<VCFImpFamilyRecord *>,
											 int>>& rs);
	void modify_00x11(
				std::vector<std::vector<VCFHeteroHomoRecord *>>& heho_recordss,
				std::vector<std::vector<VCFImpFamilyRecord *>>& other_recordss,
				const std::vector<const KnownFamily *>& families);
	std::pair<std::vector<std::vector<VCFHeteroHomoRecord *>>,
			  std::vector<std::vector<VCFImpFamilyRecord *>>>
	divide_vcf_into_record_types(const std::vector<VCFFamily *>& family_vcfs,
							const std::vector<const KnownFamily *>& families,
							const Option *option);
	std::vector<VCFFillable *> fill_vcf(
			const std::map<std::string, std::vector<VCFHeteroHomo *>>& dic_vcfs,
			const std::vector<std::vector<VCFImpFamilyRecord *>>& other_recordss,
			const std::vector<const KnownFamily *>& families, int num_threads);
	void compress_records(std::vector<VCFImpFamilyRecord *>& others);
	VCFSmall *correct_large_family_VCFs(const VCFSmall *orig_vcf,
							const std::vector<const KnownFamily *>& families,
							const Map& geno_map, const Option *option);
}

#endif
