#ifndef __VCFFILLABLERECORD
#define __VCFFILLABLERECORD

#include <vector>
#include <string>

#include "VCFImpFamily.h"
#include "ClassifyRecord.h"


//////////////////// VCFFillableRecord ////////////////////

class VCFFillableRecord : public VCFFamilyRecord {
protected:
	int			index;
	FillType	type;
	ParentComb	comb;
	
public:
	VCFFillableRecord(const STRVEC& v, const STRVEC& s,
									int i, FillType t, ParentComb c) :
						VCFFamilyRecord(v, s), index(i), type(t), comb(c) { }
	~VCFFillableRecord() { }
	
	int get_index() const { return index; }
	FillType get_type() const { return type; }
	bool is_unable() const { return type == FillType::UNABLE; }
	bool is_fillable_type() const {
		return type == FillType::UNABLE || type == FillType::IMPUTABLE;
	}
	bool is_filled_type() const {
		return type == FillType::FILLED || type == FillType::MAT
										|| type == FillType::PAT;
	}
	bool is_mat_type() const { return type == FillType::MAT; }
	bool is_pat_type() const { return type == FillType::PAT; }
	
	ParentComb get_comb() const { return comb; }
	bool is_00x11() const { return comb == ParentComb::P00x11; }
	void set_00x11() { comb = ParentComb::P00x11; }
	
	std::vector<std::pair<int, int>>	possible_phasings() const;
	
	STRVEC prog_gts() const;
	VCFFillableRecord *copy() const;
	std::tuple<int,int,int> count_gts() const;
	std::string gt_from_parent(int mat_from, int pat_from) const;
	std::string gt_from_mat(int mat_from, int c) const;
	std::string gt_from_pat(int pat_from, int c) const;
	int mat_from(int c) const;
	int pat_from(int c) const;
	
	void modify_gts();
	void modify_parents_type();
	int from_which_chrom(std::size_t i, bool is_mat) const;
	void modify_gts(const STRVEC& new_prog_gts);
	void fill_PGT();
	void set(const STRVEC& new_v, FillType new_type);
	
private:
	std::vector<std::vector<double>> make_probability_table() const;
	
	void phase();
	int find_geno_type(const std::string& type) const;
	
	std::string inverse_gt(const std::string& gt, bool inv) const;
	std::string inverse_prog_gt(const std::string& gt,
									bool inv_mat, bool inv_pat) const;
	STRVEC inverse_prog_gts(const STRVEC& prog_gts,
									bool inv_mat, bool inv_pat) const;
	void inverse_parents_gts(bool inv_mat, bool inv_pat);
	bool is_same_gts(const std::string& gt1, const std::string& gt2) const;
	bool is_near_prog_gts(const STRVEC& gts) const;
	int hash(int d) const;
	std::string decide_by_majority(const std::vector<std::string>& GTs) const;
	void swap_parents(int i, const std::string& GT);
	
public:
	static VCFRecord *integrate_records(
						const std::vector<VCFFillableRecord *>& records);
	static std::string decide_duplicated_Genotype(
							const std::vector<VCFFillableRecord *>& records,
							const std::vector<std::pair<int, int>>& positions);
	static VCFFillableRecord *convert(const VCFImpFamilyRecord *record);
	static VCFRecord *merge(const std::vector<VCFFillableRecord *>& records,
														const STRVEC& samples);
	static void integrate_each_sample(
							const std::vector<VCFFillableRecord *>& records,
							const std::vector<std::pair<int, int>>& positions);
	static bool is_all_same_GT(const std::vector<VCFFillableRecord *>& records,
						const std::vector<std::pair<int, int>>& pos_samples);
	static VCFRecord *integrate(const std::vector<VCFFillableRecord *>& records,
			const std::vector<std::string>& samples,
			const std::vector<std::vector<std::pair<int, int>>>& pos_samples);
	static int from_which_chrom(const VCFFillableRecord *record,
										std::size_t i, bool is_mat);
	static int from_which_chrom_mat(const VCFFillableRecord *record,
														std::size_t i) {
		return from_which_chrom(record, i, true);
	}
	static int from_which_chrom_pat(const VCFFillableRecord *record,
														std::size_t i) {
		return from_which_chrom(record, i, false);
	}
};
#endif
