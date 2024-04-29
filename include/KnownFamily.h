#ifndef __KNOWNFAMILY__
#define __KNOWNFAMILY__

// KnownFamily.h
#include "Pedigree.h"

// Remember if the parent is in samples
class KnownFamily : public Family {
	const bool	mat_known;
	const bool	pat_known;
	
public:
	// Use progenies after copying
	KnownFamily(const std::string& mat, const std::string& pat,
						bool mat_k, bool pat_k,
						const std::vector<const Progeny *>& progs) :
				Family(mat, pat, progs), mat_known(mat_k), pat_known(pat_k) { }
	
	~KnownFamily() { }
	
	bool is_mat_known() const { return mat_known; }
	bool is_pat_known() const { return pat_known; }
	
	bool is_any_known() const {
		return mat_known || pat_known;
	}
	
	bool is_one_unknown() const {
		return (!mat_known && pat_known) || (mat_known && !pat_known);
	}
	
	std::vector<std::string> known_parents() const {
		std::vector<std::string>	parents;
		if(this->mat_known)
			parents.push_back(this->get_mat());
		if(this->pat_known)
			parents.push_back(this->get_pat());
		return parents;
	}
	
	KnownFamily *create(const std::vector<const Progeny *>& progs) const {
		return new KnownFamily(mat, pat, mat_known, pat_known, progs);
	}
};
#endif
