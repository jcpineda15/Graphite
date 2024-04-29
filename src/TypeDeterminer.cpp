#include <set>
#include <queue>
#include <algorithm>
#include <cmath>
#include "../include/TypeDeterminer.h"

using namespace std;


//////////////////// TypeDeterminer ////////////////////

TypeDeterminer::TypeDeterminer(size_t n, double alpha_) : N(n), alpha(alpha_) {
	this->make_memo00();
	this->make_memo01();
	this->make_memo02();
	this->make_memo11();
	this->make_memo12();
	this->make_memo22();
	this->sort();
}

// cumulative probability from the center
vector<pair<int, double>> TypeDeterminer::binomial(size_t M) const {
	// Large M results in poor accuracy
	// It is better to calculate from the center
	const double	p = pow(0.5, M);
	vector<double>	ps(1, p);
	for(size_t n = 1; n <= M; ++n)
		ps.push_back(ps.back() * (M - n + 1) / n);
	
	double	total_p = 0.0;
	vector<pair<int, double>>	ps2;
	if(M % 2 == 0) {
		ps2.push_back(make_pair(M/2, total_p));
		total_p += ps[M/2];
	}
	for(int n1 = (M-1)/2; n1 >= 0; --n1) {
		const int	n2 = M - n1;
		ps2.push_back(make_pair(n1, total_p));
		ps2.push_back(make_pair(n2, total_p));
		total_p += ps[n1] * 2;
		if(total_p >= 1.0 - this->alpha)
			break;
	}
	return ps2;
}

void TypeDeterminer::insert(const State& s, ParentComb c, double p) {
	memo[s].push_back(pair<double, ParentComb>(p, c));
}

void TypeDeterminer::insert(size_t n0, size_t n1, size_t n2,
												ParentComb c, double p) {
	insert(State(n0, n1, n2), c, p);
}

void TypeDeterminer::make_memo00() {
	// Patterns that should all have the same Genotype
	// are up to 20% for all other Genotypes
	for(size_t num_NA = 0; num_NA <= N/5; ++num_NA) {
		for(size_t n1 = 0; n1 <= N/5 - num_NA; ++n1) {
			for(size_t n2 = 0; n2 <= N/5 - num_NA - n1; ++n2)
				insert(N - n1 - n2 - num_NA, n1, n2, ParentComb::P00x00, 0.0);
		}
	}
}

void TypeDeterminer::make_memo01() {
	for(size_t num_NA = 0; num_NA <= N/5; ++num_NA) {
		const size_t	L = min(N/10, N/5 - num_NA);
		for(size_t n2 = 0; n2 <= L; ++n2) {
			const int	M = N - num_NA - n2;
			const auto	ps = this->binomial(M);
			for(auto p = ps.begin(); p != ps.end(); ++p) {
				const int	n0 = p->first;
				insert(n0, M-n0, n2, ParentComb::P00x01, p->second);
			}
		}
	}
}

void TypeDeterminer::make_memo02() {
	for(size_t num_NA = 0; num_NA <= N/5; ++num_NA) {
		for(size_t n0 = 0; n0 <= N/5 - num_NA; ++n0) {
			for(size_t n2 = 0; n2 <= N/5 - num_NA - n0; ++n2)
				insert(n0, N-n0-n2-num_NA, n2, ParentComb::P00x11, 0.0);
		}
	}
}

double TypeDeterminer::genotype_probability(size_t n0,
											size_t n1, size_t n2) const {
	// If n is large, we may have to do logarithmic calculations
	const size_t	n = n0 + n1 + n2;
	const double	p0 = 0.25;
	const double	p1 = 0.5;
	const double	p2 = 0.25;
	double	p = 1.0;
	for(size_t i = 1; i <= n0; ++i)
		p *= p0 * (n - i + 1) / i;
	for(size_t i = 1; i <= n1; ++i)
		p *= p1 * (n - n0 - i + 1) / i;
	for(size_t i = 1; i <= n2; ++i)
		p *= p2 * (n2 - i + 1) / i;
	return p;
}

TypeDeterminer::PQState TypeDeterminer::initialize_state(size_t M) const {
	const size_t	n1 = M / 2;
	const size_t	n0 = (M - n1 + 1) / 2;	// let n0 >= n2
	const size_t	n2 = M - n0 - n1;
	return create_pqstate(n0, n1, n2);
}

TypeDeterminer::PQState TypeDeterminer::create_pqstate(
									size_t n0, size_t n1, size_t n2) const {
	const double	p = genotype_probability(n0, n1, n2);
	return PQState(p, n0, n1, n2);
}

TypeDeterminer::PQState TypeDeterminer::create_pqstate(const State& s) const {
	const size_t	n0 = get<0>(s);
	const size_t	n1 = get<1>(s);
	const size_t	n2 = get<2>(s);
	return create_pqstate(n0, n1, n2);
}

TypeDeterminer::State TypeDeterminer::create_state(const PQState& s) const {
	const size_t	n0 = get<1>(s);
	const size_t	n1 = get<2>(s);
	const size_t	n2 = get<3>(s);
	return State(n0, n1, n2);
}

TypeDeterminer::State TypeDeterminer::reverse_state(const State& s) const {
	const size_t	n0 = get<0>(s);
	const size_t	n1 = get<1>(s);
	const size_t	n2 = get<2>(s);
	return State(n2, n1, n0);
}

vector<TypeDeterminer::PQState> TypeDeterminer::neighbor_states(
												const PQState& s0) const {
	vector<PQState>	neighbors;
	const size_t	n1 = get<1>(s0);
	const size_t	n2 = get<2>(s0);
	const size_t	n3 = get<3>(s0);
	if(n1 > 0) {
		if(n1-1 >= n3)
			neighbors.push_back(create_pqstate(n1-1, n2+1, n3));
		if(n1-1 >= n3+1)
			neighbors.push_back(create_pqstate(n1-1, n2, n3+1));
	}
	if(n2 > 0) {
		neighbors.push_back(create_pqstate(n1+1, n2-1, n3));
		if(n1 >= n3+1)
			neighbors.push_back(create_pqstate(n1, n2-1, n3+1));
	}
	if(n3 > 0) {
		neighbors.push_back(create_pqstate(n1+1, n2, n3-1));
		neighbors.push_back(create_pqstate(n1, n2+1, n3-1));
	}
	return neighbors;
}

// 0/1 x 0/1 is difficult
void TypeDeterminer::make_memo11() {
	// sum up from the most probable state
	// until the accumulation exceeds 1-α(alpha)
	for(size_t num_NA = 0; num_NA <= N/5; ++num_NA) {
		const size_t	M = N - num_NA;
		const PQState	pqs0 = initialize_state(M);
		priority_queue<PQState>	pq;
		pq.push(pqs0);
		std::set<PQState>	visited;
		visited.insert(pqs0);
		double	total_p = 0.0;
		while(total_p < 1.0 - alpha) {
			const PQState	pqs = pq.top();
			pq.pop();
			if(get<1>(pqs) == get<3>(pqs)) {
				insert(create_state(pqs), ParentComb::P01x01, total_p);
				total_p += get<0>(pqs);
			}
			else {
				// so n0 > n2, sum up probability of n0 < n2
				const State	s1 = create_state(pqs);
				insert(s1, ParentComb::P01x01, total_p);
				insert(reverse_state(s1), ParentComb::P01x01, total_p);
				total_p += get<0>(pqs) * 2;
			}
			
			const vector<PQState>	neighbors = neighbor_states(pqs);
			for(auto p = neighbors.begin(); p != neighbors.end(); ++p) {
				const PQState	pqs1 = *p;
				if(visited.find(pqs1) == visited.end()) {
					pq.push(pqs1);
					visited.insert(pqs1);
				}
			}
		}
	}
}

void TypeDeterminer::make_memo12() {
	for(size_t num_NA = 0; num_NA <= N/5; ++num_NA) {
		const size_t	L = min(N/10, N/5 - num_NA);
		for(size_t n0 = 0; n0 <= L; ++n0) {
			const size_t	M = N - num_NA - n0;
			const auto	ps = this->binomial(M);
			for(auto p = ps.begin(); p != ps.end(); ++p) {
				const size_t	n1 = p->first;
				insert(n0, n1, M-n1, ParentComb::P01x11, p->second);
			}
		}
	}
}

void TypeDeterminer::make_memo22() {
	// same as make_memo00
	for(size_t num_NA = 0; num_NA <= N/5; ++num_NA) {
		for(size_t n0 = 0; n0 <= N/5 - num_NA; ++n0) {
			for(size_t n1 = 0; n1 <= N/5 - num_NA - n0; ++n1)
				insert(n0, n1, N - n0 - n1 - num_NA, ParentComb::P11x11, 0.0);
		}
	}
}

void TypeDeterminer::sort() {
	for(auto p = memo.begin(); p != memo.end(); ++p) {
		std::sort(p->second.begin(), p->second.end(),
							std::greater<pair<double, ParentComb>>());
	}
}

vector<pair<double, ParentComb>> TypeDeterminer::determine(
												const State& counter) const {
	const auto	p = memo.find(counter);
	if(p != memo.end())
		return p->second;
	else
		return vector<pair<double, ParentComb>>();
}

vector<int> TypeDeterminer::int_gt_pairs(int p) {
	vector<int>	v;
	for(int i = 0; i < 6; ++i) {
		if((p & (1 << i)) != 0)
			v.push_back(i);
	}
	return v;
}

pair<int,int> TypeDeterminer::int_gt_pair(ParentComb comb) {
	const int	p = static_cast<int>(comb);
	if(p == 0)
		return pair<int,int>(0, 0);
	else if(p < 3)
		return pair<int,int>(p - 1, 1);
	else
		return pair<int,int>(p - 3, 2);
}
