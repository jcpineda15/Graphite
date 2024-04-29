#include <iostream>
#include <stdexcept>
#include "../include/option.h"
#include "../include/common.h"

using namespace std;


//////////////////// Option ////////////////////

bool Option::is_efficient_chrom(int i) const {
	if(chroms.empty())
		return true;
	else
		return std::find(chroms.begin(), chroms.end(), i) != chroms.end();
}

void Option::print_info() const {
	// required
	cerr << "input VCF : " << path_vcf << endl;
	cerr << "pedigree : " << path_ped << endl;
	cerr << "output VCF : " << path_out << endl;
	
	// optional
	cerr << "number of threads : " << num_threads << endl;
	cerr << "number of progenies for large family : " << lower_progs << endl;
	
	if(!imputes_isolated_samples) {
		cerr << "isolate samples will not be imputed." << endl;
		if(outputs_unimputed_samples)
			cerr << "but, outputs these samples" << endl;
	}
}

string Option::flag_value(const string& s, int argc, char **argv) {
	for(size_t i = 1; i < (size_t)(argc - 1); ++i) {
		if(argv[i] == s)
			return argv[i+1];
	}
	return string();
}

bool Option::exists(const string& s, int argc, char **argv) {
	for(size_t i = 3U; i < (size_t)argc; ++i) {
		if(argv[i] == s)
			return true;
	}
	return false;
}

vector<size_t> Option::parse_array(const string& f) {
	vector<size_t>	array;
	const vector<string>	v = Common::split(f, ',');
	for(auto p = v.begin(); p != v.end(); ++p) {
		const vector<string>	w = Common::split(*p, ':');
		if(w.size() == 1U) {
			array.push_back(stoi(w[0]));
		}
		else if(w.size() == 2U) {
			for(size_t i = stoi(w[0]); i < (size_t)stoi(w[1]); ++i)
				array.push_back(i);
		}
		else {
			throw std::runtime_error(f + " is invalid format.");
		}
	}
	return array;
}

vector<size_t> Option::get_families(int argc, char **argv) {
	const string	s = flag_value("-f", argc, argv);
	if(s.empty())
		return vector<size_t>();
	else
		return parse_array(s);
}

vector<size_t> Option::get_chroms(int argc, char **argv) {
	const string	s = flag_value("-c", argc, argv);
	if(s.empty())
		return vector<size_t>();
	else
		return parse_array(s);
}

int Option::get_num_threads(int argc, char **argv) {
	const string	s = flag_value("-t", argc, argv);
	if(s.empty())
		return 1;
	else
		return stoi(s);
}

size_t Option::get_lower_progenies(int argc, char ** argv) {
	const string	s = flag_value("--lower-progs", argc, argv);
	if(s.empty())
		return 10;
	
	try {
		return stoul(s);
	}
	catch(const std::invalid_argument& e) {
		throw std::runtime_error(s + " is not an integer.");
	}
}

Option *Option::create(int argc, char **argv) {
	try {
		// The top three are required arguments
		const string	path_vcf = flag_value("-i", argc, argv);
		if(path_vcf.empty())
			throw std::runtime_error("input VCF not specified.");
		
		const string	path_ped = flag_value("-p", argc, argv);
		if(path_ped.empty())
			throw std::runtime_error("pedigree file not specified.");
		
		const string	path_out = flag_value("-o", argc, argv);
		if(path_out.empty())
			throw std::runtime_error("output VCF not specified.");
		
		// Optional
		const string	ref = flag_value("--ref", argc, argv);
		const string	path_map = flag_value("-m", argc, argv);
		const vector<size_t>	families = get_families(argc, argv);
		const vector<size_t>	chroms = get_chroms(argc, argv);
		const int	num_threads = get_num_threads(argc, argv);
		const size_t	lower_progs = get_lower_progenies(argc, argv);
		const bool	only_large_families = exists("-large-only", argc, argv);
		const bool	impute_isolated = !exists("--not-impute-isolated",
																argc, argv);
		const bool	out_isolated = exists("--out-isolated", argc, argv);
		if(impute_isolated && out_isolated)
			return NULL;
		
		const bool	corrects_isolated_samples
							= !exists("--correct-isolated", argc, argv);
		
		return new Option(path_vcf, ref, path_ped, path_map, families,
							chroms, num_threads, lower_progs,
							only_large_families, impute_isolated,
							out_isolated, corrects_isolated_samples, path_out);
	}
	catch(const std::runtime_error& e) {
		cerr << "error : " << e.what() << endl;
		return NULL;
	}
	catch(const std::invalid_argument& e) {
		return NULL;
	}
};

void Option::usage(char **argv) {
	cerr << argv[0] << " -i VCF [-ref ref VCF] -p ped [-m map] "
				<< "[-t num_threads] [-f family indices] [-c chrom indices] "
				<< "[--lower-progs lower num progenies] [--large-only] "
				<< "[--not-impute-isolated [--out-isolated]] "
				<< "[--correct-isolated] "
				<< "-o out." << endl;
	cerr << "family indices: (index|first:last)[,(index|first:last)[,..]]"
																	<< endl;
	cerr << "chrom indices: same as family indices." << endl;
	cerr << "--large-only: large families only." << endl;
	cerr << "--not-impute-isolated: not impute isolated samples." << endl;
	cerr << "--out-isolated: output not imputed isolated samples." << endl;
	cerr << "--correct-isolated: "
		 << "correct wrong genotypes of isolated samples." << endl;
}
