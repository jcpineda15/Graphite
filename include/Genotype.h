#ifndef __GENOTYPE
#define __GENOTYPE

#include <string>
#include <vector>
#include <utility>


//////////////////// Genotype ////////////////////

class Genotype {
	const char	gt1;
	const char	gt2;
	const bool	phasing;
	
public:
	Genotype(const std::string& s);
	
	std::pair<char,char> gts() const;
	bool includes(char gt) const;
	bool conflicts(const Genotype& mat_gt, const Genotype& pat_gt,
													bool considers_phasing);
	
	static int get_int_gt(const std::string& s);
	static std::string int_to_gt(int n);
	static bool is_NA(const std::string& gt) {
		return gt.c_str()[0] == '.' || gt.c_str()[2] == '.';
	}
	static bool is_valid(const std::string& gt, int mat_gt, int pat_gt);
	static std::string possible_gts(int gt);
	static int sum_gt(const std::string& gt);
	static bool is_all_NA(const std::vector<std::string>& GTs);
};
#endif
