#include "../include/VCFImputable.h"

using namespace std;

vector<int> VCFImputable::get_int_gts(size_t sample_index) const {
	vector<int>	int_gts;
	for(size_t i = 0; i < this->size(); ++i)
		int_gts.push_back(get_record(i)->get_int_gt(sample_index));
	return int_gts;
}

Haplotype VCFImputable::clip_haplotype(size_t sample_index, int side) const {
	const auto	hap = VCFSmallBase::clip_raw_haplotype(sample_index, side);
	return Haplotype(hap, sample_index, side);
}

HaplotypePair VCFImputable::impute_cM_each_sample(HaplotypePair prev_hap,
												size_t sample_index, bool exec,
												bool modify_genotypes) {
	const vector<int>	int_gts = get_int_gts(sample_index);
	const vector<Haplotype>	haps_mat = collect_haplotypes_mat(sample_index);
	const vector<Haplotype>	haps_pat = collect_haplotypes_pat(sample_index);
	// select a combination with seed
	const int	seed = get_record(0)->pos();
	const HaplotypePair	hap = Haplotype::impute(int_gts, haps_mat, haps_pat,
																prev_hap, seed);
	if(exec)	// actually impute?
		this->set_haplotype(hap, sample_index, modify_genotypes);
	return hap;
}

// phasing, but not correct
void VCFImputable::set_GT_unmodify(size_t i, size_t sample_id,
												const string& new_GT) {
	auto	*record = get_record(i);
	const auto	GT = record->get_GT(sample_id);
	if(record->is_NA(sample_id))
		record->set_GT(sample_id, new_GT);
	else if(GT == "0/0")
		record->set_GT(sample_id, "0|0");
	else if(GT == "1/1")
		record->set_GT(sample_id, "1|1");
	else if(stoi(new_GT.substr(0U, 1U)) + stoi(new_GT.substr(2U, 1U)) == 1)
		record->set_GT(sample_id, new_GT);
	else {
		// decide randomly whether 0|1 or 1|0
		const ll	b = 48271 * (record->pos() + sample_id) / 256 & 1;
		if(b == 0)
			record->set_GT(sample_id, "0|1");
		else
			record->set_GT(sample_id, "1|0");
	}
}

void VCFImputable::set_haplotype(HaplotypePair hap, size_t sample_index,
													bool modify_genotypes) {
	const Haplotype&	hap1 = hap.first;
	const Haplotype&	hap2 = hap.second;
	vector<string>	gts;
	for(size_t i = 0; i < hap1.hap.size(); ++i) {
		const char	gt1 = '0' + hap1.hap[i];
		const char	gt2 = '0' + hap2.hap[i];
		char	GT[] = { gt1, '|', gt2, '\0' };
		if(modify_genotypes)
			get_record(i)->set_GT(sample_index, GT);
		else
			set_GT_unmodify(i, sample_index, GT);
	}
}
