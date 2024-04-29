#include <algorithm>

#include "../include/VCFHeteroHomo.h"
#include "../include/VCFHomoHomo.h"
#include "../include/VCFHeteroHeteroLite.h"
#include "../include/VCFJunkRecord.h"
#include "../include/VCFHeteroHomoPP.h"
#include "../include/impute_prog_only.h"
#include "../include/ClassifyRecord.h"
#include "../include/option.h"

using namespace std;


//////////////////// ImputeProgOnly ////////////////////

void ImputeProgOnly::classify_record(size_t i, VCFFamilyRecord *record,
							const STRVEC& samples, const TypeDeterminer *td,
							std::vector<VCFHeteroHomoRecord *>& heho_records,
							std::vector<VCFImpFamilyRecord *>& other_records) {
	ClassifyRecord	*CR = ClassifyRecord::get_instance();
	const auto	pair1 = CR->classify(record, td, false);
	const ParentComb	pc = pair1.first;
	const WrongType	wrong_type = pair1.second;
	const auto&	v = record->get_v();
	if(TypeDeterminer::is_homohomo(pc)) {
		auto	*record_ = new VCFHomoHomoRecord(v, samples, i, wrong_type, pc);
		other_records[i] = record_->impute();
	}
	else if(TypeDeterminer::is_heterohomo(pc)) {
		heho_records[i] = new VCFHeteroHomoRecord(v, samples, i,
														wrong_type, pc);
	}
	else if(pc == ParentComb::P01x01) {
		other_records[i] = new VCFHeteroHeteroLiteRecord(v, samples, i,
																wrong_type, pc);
	}
	else {	// no candidate
		other_records[i] = new VCFJunkRecord(v, samples, i, wrong_type);
	}
}

pair<vector<VCFHeteroHomoRecord *>, vector<VCFImpFamilyRecord *>>
ImputeProgOnly::classify_records(VCFFamily *vcf,
									const STRVEC& samples, Option *option) {
	ClassifyRecord	*CR = ClassifyRecord::get_instance();
	const TypeDeterminer	*td = CR->get_TypeDeterminer(samples.size()-2,
																option->ratio);
	const size_t	M = vcf->size();
	vector<VCFHeteroHomoRecord *>	records1(M, NULL);
	vector<VCFImpFamilyRecord *>	records2(M, NULL);
	for(size_t i = 0; i < M; ++i) {
		VCFFamilyRecord	*record = vcf->get_family_record(i);
		classify_record(i, record, samples, td, records1, records2);
	}
	
	vector<VCFHeteroHomoRecord *>	heho_records;
	vector<VCFImpFamilyRecord *>	other_records;
	for(size_t i = 0; i < M; ++i) {
		auto	*r1 = records1[i];
		auto	*r2 = records2[i];
		if(r1 != NULL)
			heho_records.push_back(r1);
		else if(r2 != NULL)
			other_records.push_back(r2);
	}
	return make_pair(heho_records, other_records);
}

VCFHeteroHomoPP *ImputeProgOnly::merge_vcf(
							map<FillType, vector<VCFFillableRecord *>>& rss,
							const vector<STRVEC>& header,
							const STRVEC& samples, const Map& gmap) {
	vector<VCFFillableRecord *>	records;
	for(auto p = rss.begin(); p != rss.end(); ++p) {
		records.insert(records.end(), p->second.begin(), p->second.end());
	}
	std::sort(records.begin(), records.end(), 
				[](const VCFFillableRecord *lh, const VCFFillableRecord *rh)
				{ return lh->pos() < rh->pos(); });
	return new VCFHeteroHomoPP(header, samples, records, gmap);
}

VCFHeteroHomoPP *ImputeProgOnly::impute_prog_vcf_chr(
										const VCFSmallBase *parent_vcf,
										const VCFSmallBase *prog_vcf,
										const Map& gmap, const Option *option) {
	cout << "chr: " << parent_vcf->get_record(0)->chrom()
					<< parent_vcf->size() << " records "
					<< prog_vcf->size() << " records" << endl;
	// parent_vcfは両親のみ、prog_vcfは後代のみという前提
	auto	samples = parent_vcf->get_samples();
	const auto&	prog_samples = prog_vcf->get_samples();
	samples.insert(samples.end(), prog_samples.begin(), prog_samples.end());
	auto	*vcf = VCFHeteroHomoPP::merge(parent_vcf, prog_vcf,
													samples, gmap, option);
	vector<VCFFamilyRecord *>	records;
	for(size_t i = 0; i < vcf->size(); ++i)
		records.push_back(vcf->get_family_record(i));
	auto	rss = VCFHeteroHomoPP::classify_records(records);
	const auto	header = parent_vcf->trim_header(samples);
	auto	*mat_vcf = new VCFHeteroHomoPP(header, samples,
											rss[FillType::MAT], gmap);
	auto	*pat_vcf = new VCFHeteroHomoPP(header, samples,
											rss[FillType::PAT], gmap);
	mat_vcf->impute();
	pat_vcf->impute();
	auto	*merged_vcf = merge_vcf(rss, header, samples, gmap);
	merged_vcf->fill();
	return merged_vcf;
}
