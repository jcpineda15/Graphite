#include <stdexcept>
#include "../include/Genotype.h"

using namespace std;


//////////////////// Genotype ////////////////////

Genotype::Genotype(const string& s) : gt1(s.c_str()[0]), gt2(s.c_str()[2]),
												phasing(s.c_str()[1] == '|') { }

pair<char,char> Genotype::gts() const {
	return pair<char,char>(gt1, gt2);
}

bool Genotype::includes(char gt) const {
	return gt == gt1 || gt == gt2;
}

int Genotype::get_int_gt(const string& s) {
	try {
		return stoi(s.substr(0, 1)) + stoi(s.substr(2, 1));
	}
	catch(std::invalid_argument& e) {
		return -1;
	}
}

string Genotype::int_to_gt(int n) {
	switch(n) {
		case 0:  return "0/0";
		case 1:  return "0/1";
		case 2:  return "1/1";
		default: return "./.";
	}
}

bool Genotype::conflicts(const Genotype& mat_gt, const Genotype& pat_gt,
												bool considers_phasing) {
	if(considers_phasing && this->phasing) {
		return mat_gt.includes(gt1) && pat_gt.includes(gt2);
	}
	else {
		if(gt1 == gt2)
			return mat_gt.includes(gt1) && pat_gt.includes(gt2);
		else
			return (mat_gt.includes(gt1) && pat_gt.includes(gt2)) ||
					(pat_gt.includes(gt1) && mat_gt.includes(gt2));
	}
}

bool Genotype::is_valid(const string& gt, int mat_gt, int pat_gt) {
	if(gt.length() < 3)
		return false;
	
	const string	mat_gts = possible_gts(mat_gt);
	const string	pat_gts = possible_gts(pat_gt);
	return (mat_gts.find(gt.substr(0, 1)) != string::npos &&
			pat_gts.find(gt.substr(2, 1)) != string::npos) ||
		   (mat_gts.find(gt.substr(2, 1)) != string::npos &&
			pat_gts.find(gt.substr(0, 1)) != string::npos);
}

string Genotype::possible_gts(int gt) {
	switch(gt) {
		case  0: return "0";
		case  3: return "1";
		default: return "01";
	}
}

int Genotype::sum_gt(const string& gt) {
	return (int)((gt.c_str()[0] - '0') + (gt.c_str()[2] - '0'));
}

bool Genotype::is_all_NA(const vector<string>& GTs) {
	for(auto p = GTs.begin(); p != GTs.end(); ++p) {
		if(*p != "./.")
			return false;
	}
	return true;
}
