#ifndef __BAUMWELCHWITHFIXEDTS
#define __BAUMWELCHWITHFIXEDTS

#include <vector>
#include <string>
#include <array>


//////////////////// BaumWelch ////////////////////

class BaumWelch {
public:
	typedef std::array<double, 4>	TransitionMatrix;	// h0*2+h
	
private:
	typedef std::array<double, 6>	EmissionMatrix;		// s*2+h
	typedef std::vector<std::array<double, 2>>	Table;
	
	const std::string&	seq;
	const std::size_t	L;
	const std::vector<TransitionMatrix>	Ts;
	
	std::array<double, 2>	pi;
	EmissionMatrix			E;
	
public:
	BaumWelch(const std::string& seq_,
				const std::vector<TransitionMatrix>& Ts_);
	
	EmissionMatrix reverse_probs(const EmissionMatrix& M) const;
	std::string Viterbi(const EmissionMatrix& A) const;
	
private:
	std::size_t state(std::size_t t) const { return to_int(seq.c_str()[t]); }
	std::array<double, 2> initialize_pi();
	TransitionMatrix initialize_transition_matrix();
	EmissionMatrix initialize_emission_matrix();
	bool is_valid_Emission_Matrix(const EmissionMatrix& E) const;
	Table compute_alpha(const EmissionMatrix& A) const;
	Table compute_beta(const EmissionMatrix& A) const;
	Table compute_gamma(const Table& a, const Table& b) const;
	Table compute_parameters(const EmissionMatrix& A);
	void update_probabilities(std::array<double, 2>& pi1, EmissionMatrix& E1);
	bool is_near_parameters(const EmissionMatrix& E1,
							const std::array<double, 2>& pi1,
							const EmissionMatrix& E2,
							const std::array<double, 2>& pi2) const;
	bool converge_pi_E();
	bool converge_pi();
	std::array<double, 2> update_pi();
	EmissionMatrix reverse_emission_matrix() const;

public:
	static std::string impute(const std::string& seq,
								const std::vector<TransitionMatrix>& Ts);
	
	static std::vector<TransitionMatrix>
			modify_log(const std::vector<TransitionMatrix>& Ms);
	
private:
	static std::size_t to_int(char c) { return c == 'N' ? 2 : (size_t)(c-'0'); }
	static char to_char(std::size_t s) { return s < 2 ? '0' + s : 'N'; }
};
#endif
