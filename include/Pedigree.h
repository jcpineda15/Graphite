#ifndef __PEDIGREE
#define __PEDIGREE

#include <vector>
#include <string>
#include <set>
#include <utility>

class VCFBase;


//////////////////// Progeny ////////////////////

class Progeny {
	const std::string	family;
	const std::string	name;
	const std::string	mat;
	const std::string	pat;
	
public:
	Progeny(const std::string& f, const std::string& n,
			const std::string& m, const std::string& p) :
						family(f), name(n), mat(m), pat(p) { }
	
	const Progeny *copy() const;
	const std::string& get_name() const { return name; }
	const std::string& get_mat() const { return mat; }
	const std::string& get_pat() const { return pat; }
	std::pair<std::string, std::string> parents() const {
		return std::pair<std::string, std::string>(mat, pat);
	}
	bool with_parents() const;
	bool is_all_in_samples(const std::set<std::string>& samples) const;
	
public:
	static const Progeny *create(const std::vector<std::string>& v);
};


//////////////////// Family ////////////////////

class Family {
protected:
	const std::string	mat;
	const std::string	pat;
	const std::vector<const Progeny *>	progenies;
	const std::vector<std::string>	samples;
	
public:
	Family(const std::string& m, const std::string& p,
							const std::vector<const Progeny *>& progs) :
							mat(m), pat(p), progenies(progs),
							samples(collect_samples()) { }
	virtual ~Family();
	
	const std::string& get_mat() const { return mat; }
	const std::string& get_pat() const { return pat; }
	const std::vector<const Progeny *>&
						get_progenies() const { return progenies; }
	std::size_t num_progenies() const { return this->progenies.size(); }
	std::pair<std::string, std::string> parents() const {
		return std::pair<std::string, std::string>(mat, pat);
	}
	std::vector<std::string> collect_samples() const;
	const std::vector<std::string>& get_samples() const { return samples; }
};


//////////////////// PedigreeTable ////////////////////

class PedigreeTable {
	const std::vector<const Progeny *>	table;
	
public:
	PedigreeTable(const std::vector<const Progeny *>& progs) : table(progs) { }
	~PedigreeTable();
	
	std::size_t size() const { return table.size(); }
	const PedigreeTable *filter_with_parents(
						const std::vector<std::string>& samples) const;
	std::vector<const Family *> extract_families() const;
	std::vector<const Progeny *> get_progenies(const std::string& mat,
												 const std::string& pat) const;
	const PedigreeTable *limit_samples(
							const std::vector<std::string>& samples) const;
	std::vector<const Family *> make_families(
								const std::vector<std::string>& samples) const;
	
public:
	static const PedigreeTable *read(const std::string& path);
};
#endif
