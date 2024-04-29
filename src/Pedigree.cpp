#include <cassert>
#include "../include/common.h"
#include "../include/Pedigree.h"
#include "../include/VCF.h"

using namespace std;


//////////////////// Family ////////////////////

Family::~Family() {
	for(auto p = progenies.begin(); p != progenies.end(); ++p)
		delete *p;
}

vector<string> Family::collect_samples() const {
	vector<string>	samples;
	samples.push_back(this->mat);
	samples.push_back(this->pat);
	for(auto p = this->progenies.begin(); p != this->progenies.end(); ++p)
		samples.push_back((*p)->get_name());
	return samples;
}


//////////////////// Progeny ////////////////////

const Progeny *Progeny::copy() const {
	return new Progeny(this->family, this->name, this->mat, this->pat);
}

bool Progeny::with_parents() const {
	return mat != "0" && pat != "0";
}

bool Progeny::is_all_in_samples(const set<string>& samples) const {
	return samples.find(name) != samples.end() &&
		   samples.find(mat) != samples.end() &&
		   samples.find(pat) != samples.end();
}

const Progeny *Progeny::create(const vector<string>& v) {
	assert(v.size() >= 4U);
	return new Progeny(v[0], v[1], v[2], v[3]);
}


//////////////////// PedigreeTable ////////////////////

PedigreeTable::~PedigreeTable() {
	for(auto p = table.begin(); p != table.end(); ++p)
		delete *p;
}

const PedigreeTable *PedigreeTable::filter_with_parents(
									const vector<string>& samples) const {
	vector<const Progeny *>	progs;
	set<string>	set_samples(samples.begin(), samples.end());
	for(auto p = table.begin(); p != table.end(); ++p) {
		if((*p)->with_parents() && (*p)->is_all_in_samples(set_samples))
			progs.push_back((*p)->copy());
	}
	// not destruct Progeny in destructor
	return new PedigreeTable(progs);
}

vector<const Progeny *> PedigreeTable::get_progenies(const string& mat,
													  const string& pat) const {
	vector<const Progeny *>	children;
	set<string>	set_children;	// for avoiding duplication
	pair<string,string>	parents(mat, pat);
	for(auto p = table.begin(); p != table.end(); ++p) {
		if((*p)->parents() == parents) {
			if(set_children.find((*p)->get_name()) == set_children.end()) {
				children.push_back((*p)->copy());
				set_children.insert((*p)->get_name());
			}
		}
	}
	return children;
}

vector<const Family *> PedigreeTable::extract_families() const {
	// collect parents
	set<pair<string,string>>	s;
	for(auto p = this->table.begin(); p != this->table.end(); ++p)
		s.insert((*p)->parents());
	
	vector<const Family *>	families;
	for(auto p = s.begin(); p != s.end(); ++p) {
		const string&	mat = p->first;
		const string&	pat = p->second;
		const auto	progenies = this->get_progenies(mat, pat);
		const Family	*family = new Family(mat, pat, progenies);
		families.push_back(family);
	}
	return families;
}

vector<const Family *> PedigreeTable::make_families(
										const vector<string>& samples) const {
	map<pair<string, string>, vector<const Progeny *>>	progs;
	set<string>	used;
	for(auto p = this->table.begin(); p != this->table.end(); ++p) {
		const Progeny	*prog = *p;
		if(used.find(prog->get_name()) == used.end()) {
			progs[prog->parents()].push_back(prog->copy());
			used.insert(prog->get_name());
		}
	}
	
	vector<const Family *>	families;
	const set<string>	set_samples(samples.begin(), samples.end());
	for(auto p = progs.begin(); p != progs.end(); ++p) {
		const auto	parents = p->first;
		const auto	progenies = p->second;
		// Parent not found in samples of VCF is treated as unknown
		const string	mat = parents.first;
		const string	pat = parents.second;
		const Family	*family = new Family(mat, pat, progenies);
		families.push_back(family);
	}
	std::sort(families.begin(), families.end(),
				[](const Family *lh, const Family *rh)
				{ return lh->parents() < rh->parents(); });
	return families;
}

const PedigreeTable *PedigreeTable::limit_samples(
								const vector<string>& samples) const {
	set<string>	set_samples(samples.begin(), samples.end());
	vector<const Progeny *>	progs;
	for(auto p = table.begin(); p != table.end(); ++p) {
		if(set_samples.find((*p)->get_name()) != set_samples.end())
			progs.push_back((*p)->copy());
	}
	return new PedigreeTable(progs);
}

const PedigreeTable *PedigreeTable::read(const string& path) {
	const auto	table = Common::read_csv(path, ' ');
	vector<const Progeny *>	progs;
	for(auto p = table.begin(); p != table.end(); ++p)
		progs.push_back(Progeny::create(*p));
	return new PedigreeTable(progs);
}
