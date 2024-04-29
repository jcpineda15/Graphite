#ifndef __VCFIMPFAMILY
#define __VCFIMPFAMILY

// Abstract VCFFamily class that can be imputed

#include <vector>

#include "VCFFamily.h"
#include "ClassifyRecord.h"


//////////////////// enum ////////////////////

enum class FillType {
	MAT, PAT, FILLED, IMPUTABLE, UNABLE
};


//////////////////// VCFImpFamilyRecord ////////////////////

class VCFImpFamilyRecord : public VCFFamilyRecord {
public:
	typedef std::pair<VCFImpFamilyRecord *, int>	RecordWithPos;
	
protected:
	const int	index;
	WrongType	wrong_type;
	const ParentComb	comb;
	
public:
	VCFImpFamilyRecord(const STRVEC& v, const STRVEC& samples,
							int i, WrongType type, ParentComb c) :
									VCFFamilyRecord(v, samples), index(i),
									wrong_type(type), comb(c) { }
	virtual ~VCFImpFamilyRecord() { }
	
	int get_index() const { return index; }
	ParentComb get_comb() const { return comb; }
	
	bool is_right() const {
		return this->wrong_type == WrongType::RIGHT;
	}
	bool is_fixed() const {
		return (this->comb == ParentComb::P00x00 ||
					this->comb == ParentComb::P11x11) ||
				((this->comb == ParentComb::P00x01 ||
					this->comb == ParentComb::P01x11) && is_right());
	}
	bool is_fillable() const {
		return wrong_type != WrongType::MIX &&
				wrong_type != WrongType::UNSPECIFIED;
	}
	
	bool is_00x11() const { return this->comb == ParentComb::P00x11; }
	
	virtual bool is_homohomo() const = 0;
	virtual bool is_imputable() const = 0;
	virtual FillType get_fill_type() const = 0;
	
	void enable_modification() {
		this->wrong_type = WrongType::MODIFIABLE;
	}
	void set_00x11_parents(int i, int gt);
	
public:
	static int get_records_type(const std::vector<RecordWithPos>& records);
	static std::size_t find_fixed_index(
		const std::vector<std::pair<int, std::vector<RecordWithPos>>>& items);
	static std::pair<int, std::vector<RecordWithPos>>
	which_is_fixed(const std::vector<RecordWithPos>& v);
	static void modify_00x11(const std::vector<VCFImpFamilyRecord *>& records);
};
#endif
