#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "../include/log.h"
#include "../include/Baum_Welch_with_fixed_Ts.h"

using namespace std;


BaumWelch::BaumWelch(const string& seq_, const vector<TransitionMatrix>& Ts_) :
											seq(seq_),
											L(seq_.size()),
											Ts(modify_log(Ts_)),
											pi(initialize_pi()),
											E(initialize_emission_matrix()) {
	if(!converge_pi_E())
		converge_pi();
}

vector<BaumWelch::TransitionMatrix> BaumWelch::modify_log(
										const vector<TransitionMatrix>& Ms_) {
	vector<TransitionMatrix>	Ms;
	for(auto p = Ms_.begin(); p != Ms_.end(); ++p) {
		TransitionMatrix	M;
		for(size_t i = 0; i < p->size(); ++i)
			M[i] = Log::modified_log((*p)[i]);
		Ms.push_back(M);
	}
	return Ms;
}

std::array<double, 2> BaumWelch::initialize_pi() {
	// dividing equally
	const double	prob = -std::log(2.0);
	array<double, 2>	probs { prob, prob };
	return probs;
}

BaumWelch::TransitionMatrix BaumWelch::initialize_transition_matrix() {
	TransitionMatrix	T;
	for(size_t h1 = 0; h1 < 2; ++h1) {
		for(size_t h2 = 0; h2 < 2; ++h2) {
			T[h1*2+h2] = std::log(h1 == h2 ? 0.9 : 0.1);
		}
	}
	return T;
}

BaumWelch::EmissionMatrix BaumWelch::initialize_emission_matrix() {
	// probability of hidden state -> state
	EmissionMatrix	E;
	const size_t	N = 3;
	for(size_t h = 0; h < 2; ++h) {
		for(size_t s = 0; s < 3; ++s) {
			E[s*2+h] = std::log(h == s ? 0.9 : 0.1 / (N - 1));
		}
	}
	return E;
}

BaumWelch::Table BaumWelch::compute_alpha(const EmissionMatrix& A) const {
	Table	a(L);	// probability of s
	for(size_t h = 0; h < 2; ++h) {
		a[0][h] = pi[h] + A[state(0)*2+h];
	}
	
	for(size_t t = 1U; t < seq.size(); ++t) {
		const TransitionMatrix&	T = Ts[t-1];
		const size_t	s = state(t);
		for(size_t h = 0; h < 2; ++h) {
			array<double, 2>	v;
			for(size_t h_prev = 0; h_prev < 2; ++h_prev)
				v[h_prev] = a[t-1][h_prev] + T[h_prev*2+h] + A[s*2+h];
			a[t][h] = Log::add(v[0], v[1]);
		}
	}
	return a;
}

BaumWelch::Table BaumWelch::compute_beta(const EmissionMatrix& A) const {
	Table	b(L);
	b[L-1] = array<double, 2> { 0.0, 0.0 };
	
	for(int t = (int)L - 2; t >= 0; --t) {
		const TransitionMatrix&	T = Ts[t];
		const size_t	s = state(t + 1);
		for(size_t h = 0; h < 2; ++h) {
			array<double, 2>	v;
			for(size_t h_next = 0; h_next < 2; ++h_next)
				v[h_next] = T[h*2+h_next] + A[s*2+h_next] + b[t+1][h_next];
			b[t][h] = Log::add(v[0], v[1]);
		}
	}
	return b;
}

BaumWelch::Table BaumWelch::compute_gamma(
						const Table& a, const Table& b) const {
	Table	gamma1(L);
	for(size_t t = 0; t < L; ++t) {
		for(size_t h = 0; h < 2; ++h)
			gamma1[t][h] = a[t][h] + b[t][h];
	}
	
	// probability that seq is observed
	const double	prob_seq = Log::add(gamma1[0][0], gamma1[0][1]);
	
	Table	gamma(L);
	for(size_t t = 0; t < L; ++t) {
		for(size_t h = 0; h < 2; ++h)
			gamma[t][h] = gamma1[t][h] - prob_seq;
	}
	
	return gamma;
}

BaumWelch::Table BaumWelch::compute_parameters(const EmissionMatrix& A) {
	// probability of s
	const Table	a = compute_alpha(A);
	const Table	b = compute_beta(A);
	return compute_gamma(a, b);
}

void BaumWelch::update_probabilities(array<double, 2>& pi1,
										EmissionMatrix& E1) {
	const EmissionMatrix	A = reverse_probs(E);
	Table	gamma = compute_parameters(A);
	pi1 = gamma[0];
	
	for(size_t hj = 0; hj < 2; ++hj) {
		vector<double>	v(L);
		for(size_t t = 0; t < L; ++t)
			v[t] = gamma[t][hj];
		const double	cum = Log::sum(v);
		
		array<vector<double>, 3>	ws;
		for(size_t t = 0; t < L; ++t) {
			const size_t	s = state(t);
			ws[s].push_back(gamma[t][hj]);
		}
		
		for(size_t s = 0; s < 3; ++s) {
			if(!ws[s].empty())
				E1[s*2+hj] = Log::sum(ws[s]) - cum;
		}
	}
}

bool BaumWelch::is_valid_Emission_Matrix(const EmissionMatrix& E) const {
	for(size_t h = 0; h < 2; ++h) {
		if(E[h*2+h] <= std::log(0.5))
			return false;
	}
	return true;
}

bool BaumWelch::is_near_parameters(const EmissionMatrix& E1,
								   const array<double, 2>& pi1,
								   const EmissionMatrix& E2,
								   const array<double, 2>& pi2) const {
	vector<double>	v1(6);
	for(size_t i = 0; i < 6; ++i) {
		if(E1[i] == 0.0 || E2[i] == 0.0)
			v1[i] = Log::LOGZERO;
		else
			v1[i] = Log::diff(E1[i], E2[i]);
	}
	const double	d2 = Log::sum(v1);
	
	const double	d3 = Log::add(Log::diff(pi1[0], pi2[0]),
								  Log::diff(pi1[1], pi2[1]));
	
	return Log::add(d2, d3) < std::log(1e-7);
}

bool BaumWelch::converge_pi_E() {
	for(int i = 0; i < 100; ++i) {
		try {
			array<double, 2>	pi1;
			// if seq has not 'N', some elements of E are not updated
			EmissionMatrix		E1 = {0.0};
			update_probabilities(pi1, E1);
			if(is_near_parameters(E, pi, E1, pi1))
				break;
			pi = pi1;
			E =  E1;
		}
		catch(runtime_error& e) {
			break;
		}
	}
	
	return is_valid_Emission_Matrix(E);
}

BaumWelch::EmissionMatrix BaumWelch::reverse_probs(
										const EmissionMatrix& M) const {
#if 1
	EmissionMatrix	rev_M;	// index: s*2+h
	for(size_t s = 0; s < 3; ++s) {
		double	sum_observed = Log::add(M[s*2], M[s*2+1]);
		for(size_t h = 0; h < 2; ++h) {
			rev_M[s*2+h] = M[s*2+h] - sum_observed;
		}
	}
#else
	array<double, 3>	sum_observed;
	for(size_t s0 = 0; s0 < 3; ++s0) {
		sum_observed[s0] = Log::add(M[s0*2], M[s0*2+1]);
	}
	
	EmissionMatrix	rev_M;	// index: s*2+h
	for(size_t h = 0; h < 2; ++h) {
		for(size_t s = 0; s < 3; ++s)
			rev_M[s*2+h] = M[s*2+h] - sum_observed[s];
	}
#endif
	return rev_M;
}

array<double, 2> BaumWelch::update_pi() {
	const auto	A = reverse_probs(E);
	const auto	gamma = compute_parameters(A);
	const auto	pi1 = gamma.front();
	return pi1;
}

bool BaumWelch::converge_pi() {
	pi = initialize_pi();
	E = initialize_emission_matrix();
	for(int i = 0; i < 100; ++i) {
		try {
			const auto	pi1 = update_pi();
			double	log_s = Log::add(Log::diff(pi[0], pi1[0]),
									 Log::diff(pi[1], pi1[1]));
			if(log_s < std::log(1e-7))
				break;
			pi = pi1;
		}
		catch(runtime_error& e) {
			break;
		}
	}
	return true;
}

string BaumWelch::Viterbi(const EmissionMatrix& A) const {
	Table	table(L);
	for(size_t h = 0; h < 2; ++h) {
		table[0][h] = pi[h] + A[state(0)*2+h];
	}
	
	for(size_t t = 1; t < L; ++t) {
		const size_t	s = state(t);
		const auto&	T = Ts[t-1];
		for(size_t h = 0; h < 2; ++h) {
			double max_value = Log::LOGZERO;
			for(size_t h0 = 0; h0 < 2; ++h0) {
				const double	v = table[t-1][h0] + T[h0*2+h] + A[s*2+h];
				max_value = std::max(max_value, v);
			}
			table[t][h] = max_value;
		}
	}
	
	for(size_t t = 1; t < L; ++t) {
		const size_t	s = state(t);
		const auto&		T = Ts[t-1];
		for(size_t h = 0; h < 2; ++h) {
			double	max_value = Log::LOGZERO;
			for(size_t h0 = 0; h0 < 2; ++h0) {
				const double	v = table[t-1][h0] + T[h0*2+h] + A[s*2+h];
				max_value = std::max(max_value, v);
			}
			table[t][h] = max_value;
		}
	}
	
	// backtrack
	string	new_hidden_seq = "";
	const auto&	v = table.back();
	size_t	h = v[0] >= v[1] ? 0 : 1;
	new_hidden_seq += to_char(h);
	
	for(size_t t = L-1; t > 0; --t) {
		const size_t	s = state(t);
		const TransitionMatrix&	T = Ts[t-1];
		const double	v1 = table[t-1][0] + T[0*2+h] + A[s*2+h];
		const double	v2 = table[t-1][1] + T[1*2+h] + A[s*2+h];
		h = v1 >= v2 ? 0 : 1;
		new_hidden_seq += to_char(h);
	}
	
	std::reverse(new_hidden_seq.begin(), new_hidden_seq.end());
	return new_hidden_seq;
}

BaumWelch::EmissionMatrix BaumWelch::reverse_emission_matrix() const {
	return this->reverse_probs(E);
}

string BaumWelch::impute(const string& seq,
							const vector<TransitionMatrix>& Ts) {
	BaumWelch	BW(seq, Ts);
	const auto	A = BW.reverse_emission_matrix();
	return BW.Viterbi(A);
}
