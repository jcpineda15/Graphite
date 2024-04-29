// common.h
#ifndef __COMMON
#define __COMMON

#include <vector>
#include <string>
#include <algorithm>

namespace Common {
	std::string strip(const std::string& s);
	bool empty_line(const std::string& s);
	std::vector<std::string> split(const std::string& s, char delim);
	std::string join(const std::vector<std::string>& v, char delim);
	void chop(char *buff);
	std::vector<std::string> merge_vector(const std::vector<std::string>& v1,
										  const std::vector<std::string>& v2);
	
	std::vector<std::vector<std::string>> read_csv(
							const std::string& path, char delim=',');
	std::vector<std::vector<std::string>> read_tsv(const std::string& path);
	void write_tsv(const std::vector<std::string>& v, std::ostream& os);
	
	// e.g. [2, 5, 7], [3, 5, 6] -> [2, 5, 7, 3, 6]
	std::vector<std::string> merge_vectors(
								const std::vector<std::string>& v1,
								const std::vector<std::string>& v2);
	bool is_all_same(const std::string& seq);
	template<typename T>
	bool is_all_same(const std::vector<T>& v) {
		for(auto p = v.begin() + 1; p != v.end(); ++p) {
			if(*p != v.front())
				return false;
		}
		return true;
	}
	
	template<typename T>
	std::vector<T> unique_vector(const std::vector<T>& v) {
		std::vector<T>	w(1U, v.front());
		for(auto p = v.begin() + 1; p != v.end(); ++p) {
			if(std::find(w.begin(), w.end(), *p) == w.end())
				w.push_back(*p);
		}
		return w;
	}
	
	template<typename T>
	void delete_all(const std::vector<T *>& v) {
		for(auto p = v.begin(); p != v.end(); ++p)
			delete *p;
	}
}
#endif
