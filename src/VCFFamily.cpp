#include <algorithm>
#include <cassert>
#include "../include/VCFFamily.h"

using namespace std;


//////////////////// VCFFamilyRecord ////////////////////

VCFFamilyRecord *VCFFamilyRecord::copy() const {
	return new VCFFamilyRecord(v, samples);
}

vector<int> VCFFamilyRecord::progeny_gts() const {
	vector<int>	w(this->num_progenies());
	for(size_t i = 2U; i < samples.size(); ++i)
		w[i-2] = this->get_int_gt(i);
	return w;
}

vector<int> VCFFamilyRecord::get_progeny_int_gts() const {
	vector<int>	w(this->num_progenies());
	for(size_t i = 2U; i < this->samples.size(); ++i) {
		assert(i < w.size() + 2);
		w[i-2] = get_int_gt(i);
	}
	return w;
}

tuple<int,int,int> VCFFamilyRecord::count_gts() const {
	int	counter[] = { 0, 0, 0 };
	const auto	gts = progeny_gts();
	for(auto p = gts.begin(); p != gts.end(); ++p) {
		if(0 <= *p && *p <= 2)
			counter[*p] += 1;
	}
	return make_tuple(counter[0], counter[1], counter[2]);
}

void VCFFamilyRecord::set(const STRVEC& new_v) {
	v = new_v;
	VCFRecord::set(v);
}

void VCFFamilyRecord::impute_homohomo() {
	const char	mat_gt = this->mat_gt().c_str()[0];
	const char	pat_gt = this->pat_gt().c_str()[0];
	char	chars_GT[4] = { mat_gt, '|', pat_gt, '\0' };
	const string	GT(chars_GT);
	for(size_t i = 2; i != samples.size(); ++i)
		this->set_GT(i, GT);
}


//////////////////// VCFFamily ////////////////////

VCFFamily::VCFFamily(const vector<STRVEC>& h, const STRVEC& s,
										vector<VCFFamilyRecord *> rs) :
				VCFBase(h, s), VCFSmallBase(), VCFFamilyBase(), records(rs) { }

VCFFamily::~VCFFamily() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

VCFFamily *VCFFamily::create(const VCFSmall *vcf, const STRVEC& samples) {
	const auto	columns = vcf->extract_columns(samples);
	const auto	header = vcf->trim_header(samples);
	const vector<VCFRecord *>&	orig_records = vcf->get_records();
	vector<VCFFamilyRecord *>	records;
	for(auto p = orig_records.begin(); p != orig_records.end(); ++p)
		records.push_back(VCFFamily::subset(*p, samples, columns));
	VCFFamily	*new_vcf = new VCFFamily(header, samples, records);
	return new_vcf;
}

VCFFamilyRecord *VCFFamily::subset(VCFRecord *record, const STRVEC& samples,
												const vector<size_t>& columns) {
	const STRVEC&	orig_v = record->get_v();
	STRVEC	v(orig_v.begin(), orig_v.begin() + 9);
	for(auto p = columns.begin(); p != columns.end(); ++p) {
		v.push_back(*p != string::npos ? orig_v[*p] : "./.");
	}
	return new VCFFamilyRecord(v, samples);
}

VCFFamily *VCFFamily::merge(const VCFFamily *vcf1, const VCFFamily *vcf2) {
	// Create an empty VCF to make samples common to VFC and Record
	// and use the samples
	vector<VCFFamilyRecord *>	records;
	VCFFamily	*vcf = new VCFFamily(vcf1->get_header(),
											vcf1->get_samples(), records);
	vcf1->copy_chrs(vcf);
	
	size_t	k = 0U;
	size_t	l = 0U;
	while(k < vcf1->size() && l < vcf2->size()) {
		VCFFamilyRecord	*record1 = vcf1->get_family_record(k);
		VCFFamilyRecord	*record2 = vcf2->get_family_record(l);
		POSITION	pos1 = vcf1->record_position(*record1);
		POSITION	pos2 = vcf2->record_position(*record2);
		assert(pos1 != pos2);
		if(pos1 < pos2) {
			records.push_back(record1->copy());
			k += 1;
		}
		else {
			records.push_back(record2->copy());
			l += 1;
		}
	}
	
	for(size_t i = k; i < vcf1->size(); ++i)
		records.push_back(vcf1->get_family_record(i)->copy());
	for(size_t i = l; i < vcf2->size(); ++i)
		records.push_back(vcf2->get_family_record(i)->copy());
	
	vcf->set_records(records);
	return vcf;
}

VCFFamily *VCFFamily::create_by_two_vcfs(const VCFSmallBase *vcf1,
										 const VCFSmallBase *vcf2,
										 const STRVEC& samples) {
	// assume that samples are [mat, pat, prog1, prog2, ...]
	const vector<size_t>	columns1 = vcf1->extract_columns(samples);
	const vector<size_t>	columns2 = vcf2->extract_columns(samples);
	const auto	new_header = vcf1->trim_header(samples);
	
	vector<VCFFamilyRecord *>	new_records;
	for(size_t i = 0; i < vcf1->size(); ++i) {
		const auto	*record1 = vcf1->get_record(i);
		const auto	*record2 = vcf2->get_record(i);
		STRVEC	v(record1->get_v().begin(), record1->get_v().begin() + 9);
		for(size_t i = 0; i < samples.size(); ++i) {
			if(columns1[i] != string::npos)
				v.push_back(record1->get_v()[columns1[i]]);
			else if(columns2[i] != string::npos)
				v.push_back(record2->get_v()[columns2[i]]);
			else
				v.push_back("./.");
		}
		auto	*new_record = new VCFFamilyRecord(v, samples);
		new_records.push_back(new_record);
	}
	return new VCFFamily(new_header, samples, new_records);
}

VCFFamily *VCFFamily::join(const vector<VCFFamily *>& vcfs) {
	vector<VCFFamilyRecord *>	records;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		VCFFamily	*vcf = *p;
		for(size_t i = 0U; i < vcf->size(); ++i) {
			VCFFamilyRecord	*record = vcf->get_family_record(i);
			records.push_back(record->copy());	// for memory leak
		}
	}
	VCFFamily	*vcf = vcfs.front();
	return new VCFFamily(vcf->get_header(), vcf->get_samples(), records);
}
