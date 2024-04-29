// common.cpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "../include/common.h"

using namespace std;

string Common::strip(const string& s) {
	size_t	p1 = 0U;
	for( ; p1 < s.length(); ++p1) {
		const char	c = s.c_str()[p1];
		if(c != ' ' && c != '\t')
			break;
	}
	if(p1 == s.length())
		return "";
	size_t	p2 = s.length() - 1;
	for( ; p2 > p1; --p2) {
		const char	c = s.c_str()[p2];
		if(c != ' ' && c != '\t')
			break;
	}
	
	return s.substr(p1, p2 - p1 + 1);
}

bool Common::empty_line(const string& s) {
	const char	whitespace[] = " \t\r\n";
	for(auto p = s.begin(); p != s.end(); ++p) {
		size_t q = 0U;
		for( ; q < 4U; ++q) {
			if(*p == whitespace[q])
				break;
		}
		if(q == 4U)
			return false;
	}
	return true;
}

vector<string> Common::split(const string& s, char delim) {
	if(s.empty())
		return vector<string>();
	
	vector<string>	v;
	size_t	first = 0;
	for(size_t p = 0; p < s.size(); ++p) {
		if(s.c_str()[p] == delim) {
			if(p == first)
				v.push_back("");
			else
				v.push_back(s.substr(first, p-first));
			first = p + 1;
		}
	}
	v.push_back(s.substr(first, s.length()-first));
	return v;
}

string Common::join(const vector<string>& v, char delim) {
	stringstream	ss;
	if(!v.empty())
		ss << v.front();
	for(auto p = v.begin() + 1; p != v.end(); ++p)
		ss << delim << *p;
	return ss.str();
}

void Common::chop(char *buff) {
	char	*p;
	for(p = buff; *p != '\0'; ++p)
		;
	
	--p;
	if(*p != '\n' && *p != '\r')
		return;
	*p = '\0';
	
	--p;
	if(*p != '\n' && *p != '\r')
		return;
	*p = '\0';
}

vector<string> Common::merge_vector(const vector<string>& v1,
									const vector<string>& v2) {
	vector<string>	v = v1;
	for(auto p = v2.begin(); p != v2.end(); ++p) {
		auto	q = find(v.begin(), v.end(), *p);
		if(q == v.end())
			v.push_back(*p);
	}
	return v;
}

vector<vector<string>> Common::read_csv(const string& path, char delim) {
	ifstream	ifs(path.c_str());
	if(!ifs) {
		cerr << "error : can't open " << path << "." << endl;
		exit(1);
	}
	
	vector<vector<string>>	table;
	string	line;
	while(getline(ifs, line)) {
		if(line.c_str()[line.length()-1] == '\r')
			line = line.substr(0, line.length()-1);
		auto	v = split(line, delim);
		table.push_back(v);
	}
	return table;
}

vector<vector<string>> Common::read_tsv(const string& path) {
	return read_csv(path, '\t');
}

void Common::write_tsv(const vector<string>& v, ostream& os) {
	if(!v.empty()) {
		os << v.front();
		for(auto p = v.begin() + 1; p != v.end(); ++p)
			os << '\t' << *p;
	}
	os << endl;
}

vector<string> Common::merge_vectors(const vector<string>& v1,
									 const vector<string>& v2) {
	vector<string>	v = v1;
	for(auto p = v2.begin(); p != v2.end(); ++p) {
		if(find(v.begin(), v.end(), *p) == v.end())
			v.push_back(*p);
	}
	return v;
}

bool Common::is_all_same(const string& seq) {
	for(auto p = seq.begin() + 1; p != seq.end(); ++p) {
		if(*p != seq.front())
			return false;
	}
	return true;
}
