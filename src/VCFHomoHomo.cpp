#include "../include/common.h"
#include "../include/VCFHomoHomo.h"

using namespace std;


//////////////////// VCFHomoHomoRecord ////////////////////

bool VCFHomoHomoRecord::is_imputable() const {
	return this->comb == ParentComb::P00x11 && this->get_GT(0) == "./.";
}

vector<string> VCFHomoHomoRecord::gts() const {
	if(this->comb == ParentComb::P00x00)
		return vector<string>(samples.size(), "0|0");
	else if(this->comb == ParentComb::P11x11)
		return vector<string>(samples.size(), "1|1");
	
	// 0/0 x 1/1 is difficult
	if((this->mat_int_gt() == 0 && this->pat_int_gt() != 0) ||
				(this->mat_int_gt() != 2 && this->pat_int_gt() == 2)) {
		vector<string>	gts(samples.size(), "0|1");
		gts[0] = "0|0";
		gts[1] = "1|1";
		return gts;
	}
	else if((this->mat_int_gt() == 2 && this->pat_int_gt() != 2) ||
				(this->mat_int_gt() != 0 && this->pat_int_gt() == 0)) {
		vector<string>	gts(samples.size(), "1|0");
		gts[0] = "1|1";
		gts[1] = "0|0";
		return gts;
	}
	else {
		vector<string>	gts(samples.size(), "0/1");
		gts[0] = "./.";
		gts[1] = "./.";
		return gts;
	}
}

VCFHomoHomoRecord *VCFHomoHomoRecord::impute() const {
	vector<string>	w = v;
	const vector<string>	gts1 = this->gts();
	for(size_t i = 9; i < w.size(); ++i)
		w[i] = gts1[i-9] + w[i].substr(3);
	return new VCFHomoHomoRecord(w, this->samples, this->index,
									this->wrong_type, this->comb);
}


//////////////////// VCFHomoHomo ////////////////////

VCFHomoHomo::VCFHomoHomo(const vector<STRVEC>& h, const STRVEC& s,
									vector<VCFHomoHomoRecord *> rs) :
							VCFBase(h, s), VCFFamilyBase(), records(rs) { }

vector<VCFHomoHomo *> VCFHomoHomo::impute() const {
	vector<VCFHomoHomoRecord *>	records;
	for(auto p = records.begin(); p != records.end(); ++p)
		records.push_back((*p)->impute());
	// HeteroHetero returns a vector and matches it
	return vector<VCFHomoHomo *>(1U,
						new VCFHomoHomo(this->header, this->samples, records));
}
