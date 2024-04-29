#include <set>
#include <stack>
#include <algorithm>
#include <cassert>
#include <random>
#include "../include/inverse_graph.h"

using namespace std;


//////////////////// InverseGraph ////////////////////

InverseGraph *InverseGraph::copy() const {
	auto	*g = new InverseGraph();
	for(auto p = begin(); p != end(); ++p)
		(*g)[p->first] = p->second;
	return g;
}

vector<InverseGraph::Node> InverseGraph::collect_nodes() const {
	vector<Node>	nodes;
	for(auto p = begin(); p != end(); ++p)
		nodes.push_back(p->first);
	return nodes;
}

vector<InverseGraph::Node> InverseGraph::neighbors(Node v0) const {
	auto	p = find(v0);
	if(p == end())
		return vector<Node>();
	else {
		vector<Node>	neighs;
		for(auto q = p->second.begin(); q != p->second.end(); ++q)
			neighs.push_back(get<0>(*q));
		return neighs;
	}
}

vector<InverseGraph::Node> InverseGraph::sort_nodes() const {
	vector<Node>	vs;
	for(auto p = begin(); p != end(); ++p)
		vs.push_back(p->first);
	return vs;
}

map<InverseGraph::Node, size_t> InverseGraph::calc_dic() const {
	map<Node, size_t>	dic;
	size_t	i = 0;
	for(auto p = begin(); p != end(); ++i, ++p)
		dic[p->first] = i;
	return dic;
}

vector<bool> InverseGraph::optimize_inversions() const {
	const auto	subgraphs = this->divide_into_connected();
	vector<bool>	bs(this->size(), false);
	for(auto p = subgraphs.begin(); p != subgraphs.end(); ++p) {
		const InverseGraph&	subg = *p;
		const auto	dic_bs = subg.optimize_inversions_connected();
		for(auto q = dic_bs.begin(); q != dic_bs.end(); ++q)
			bs[q->first] = q->second;
	}
	return bs;
}

vector<InverseGraph> InverseGraph::divide_into_connected() const {
	const vector<vector<Node>>	nss = divide_nodes_into_connected();
	vector<InverseGraph>	graphs(nss.size());
	for(size_t i = 0; i < nss.size(); ++i) {
		for(auto p = nss[i].begin(); p != nss[i].end(); ++p) {
			const Node	v = *p;
			auto	q = find(v);
			graphs[i][v] = q->second;
		}
	}
	return graphs;
}

bool InverseGraph::is_leaf(Node v) const {
	auto	p = find(v);
	return p != end() && p->second.size() == 1;
}

void InverseGraph::remove_directed_edge(Node v1, Node v2) {
	const auto	neighs1 = this->neighbors(v1);
	if(neighs1.size() == 1)
		this->erase(v1);
	else {
		auto&	vs2 = (*this)[v1];
		for(auto p = vs2.begin(); p != vs2.end(); ++p) {
			if(get<0>(*p) == v2) {
				vs2.erase(p, p+1);
				return;
			}
		}
		throw("error in InverseGraph::remove_directed_edge");
	}
}

void InverseGraph::remove_edge(Node v1, Node v2) {
	this->remove_directed_edge(v1, v2);
	this->remove_directed_edge(v2, v1);
}

pair<InverseGraph *, vector<InverseGraph::Node>>
					InverseGraph::remove_acyclic_nodes() const {
	auto	*graph = this->copy();
	vector<Node>	removed_nodes;
	while(true) {
		const size_t	prev_num_removed_nodes = removed_nodes.size();
		const auto	vs = graph->collect_nodes();
		for(auto p = vs.begin(); p != vs.end(); ++p) {
			const Node	v = *p;
			if(graph->is_leaf(v)) {
				const auto	neighs = graph->neighbors(v);
				const Node	v1 = neighs[0];
				graph->remove_edge(v, v1);
				removed_nodes.push_back(v);
			}
		}
		if(removed_nodes.size() == prev_num_removed_nodes)
			break;
	}
	return make_pair(graph, removed_nodes);
}

// assume connected
bool InverseGraph::is_tree() const {
	size_t	num_edges = 0;
	for(auto p = begin(); p != end(); ++p)
		num_edges += p->second.size();
	return num_edges == (size() - 1) * 2;
}

map<InverseGraph::Node, bool> InverseGraph::decide_inversions_tree() const {
	map<InverseGraph::Node, bool>	dic_bs;
	const Node	v0 = begin()->first;
	set<Node>	visited;
	visited.insert(v0);
	vector<pair<Node, bool>>	stk(1, make_pair(v0, false));
	dic_bs[v0] = false;
	while(!stk.empty()) {
		const auto	_p = stk.back();
		stk.pop_back();
		const Node	v = _p.first;
		const bool	b = _p.second;
		auto	p = find(v);
		if(p == end())
			continue;
		for(auto q = p->second.begin(); q != p->second.end(); ++q) {
			const Node	v1 = get<0>(*q);
			const int	n1 = get<1>(*q);
			const int	n2 = get<2>(*q);
			const bool	b1 = b ^ (n1 < n2);
			if(visited.find(v1) == visited.end()) {
				dic_bs[v1] = b1;
				stk.push_back(make_pair(v1, b1));
				visited.insert(v1);
			}
		}
	}
	return dic_bs;
}

map<InverseGraph::Node, bool>
				InverseGraph::optimize_inversions_connected_core() const {
	if(this->size() <= 20)
		return this->search_all();
	
	// Determine edge from those with large true/false bias
	// It is OK if there is no contradiction until all the graphs are connected.
	const auto	dic_bs = this->connect_biased_edges();
	if(!dic_bs.empty())
		return dic_bs;
	
	// decide to flip randomly if contradictory
	// I really want to use SA
	return this->search_randomly();
}

void InverseGraph::add_removed_nodes(map<InverseGraph::Node, bool>& dic_bs,
									const vector<Node>& removed_nodes) const {
	for(auto p = removed_nodes.rbegin(); p != removed_nodes.rend(); ++p) {
		const Node	v1 = *p;
		const auto	vs = find(v1)->second;
		for(auto q = vs.begin(); q != vs.end(); ++q) {
			const Node	v2 = get<0>(*q);
			const int	n1 = get<1>(*q);
			const int	n2 = get<2>(*q);
			auto	r = dic_bs.find(v2);
			if(r != dic_bs.end()) {
				const bool	b2 = r->second;
				const bool	b1 = b2 ^ (n1 < n2);
				dic_bs.insert(make_pair(v1, b1));
				break;
			}
		}
	}
}

// graph must be connected
map<InverseGraph::Node, bool>
					InverseGraph::optimize_inversions_connected() const {
	if(this->is_tree())
		return this->decide_inversions_tree();
	
	const auto	p = this->remove_acyclic_nodes();
	auto	*graph = p.first;
	auto	removed_nodes = p.second;
	auto	dic_bs = graph->optimize_inversions_connected_core();
	this->add_removed_nodes(dic_bs, removed_nodes);
	delete graph;
	return dic_bs;
}

// brute force
map<InverseGraph::Node, bool> InverseGraph::search_all() const {
	const size_t	N = this->size();
	int	min_score = 100000000;
	vector<bool>	min_invs;
	for(int i = 0; i < (1 << (N-1)); ++i) {
		// Instead of converting it to bools,
		// it's faster to keep it as an integer.
		vector<bool>	invs(N-1);
		for(size_t j = 0; j < N - 1; ++j)
			invs[j] = ((i >> j) & 1) == 1;
		
		const int	score = this->calc_match_score(invs);
		if(score < min_score) {
			min_score = score;
			min_invs = invs;
		}
	}
	
	map<Node, bool>	dic;
	const vector<Node>	vs = sort_nodes();
	dic[vs[0]] = false;
	for(size_t i = 1; i < vs.size(); ++i)
		dic[vs[i]] = min_invs[i-1];
	return dic;
}

int InverseGraph::choice_seed() const {
	int	s = 0;
	int	D = 1000000000;
	for(auto p = begin(); p != end(); ++p) {
		for(auto q = p->second.begin(); q != p->second.end(); ++q)
			s = (s + (int)(get<0>(*q))) + D;
	}
	return s;
}

std::map<InverseGraph::Node, bool> InverseGraph::search_randomly() const {
	vector<bool>	min_bs;
	int				min_score = 1000000000;
	std::mt19937 generator(this->choice_seed());
	std::bernoulli_distribution	distribution(0.5);
	for(int i = 0; i < 524288; ++i) {	// 2^19
		vector<bool>	bs(this->size() - 1);
		for(size_t k = 0; k < this->size() - 1; ++k)
			bs.push_back(distribution(generator));
		const int	cur_score = this->calc_match_score(bs);
		if(min_bs.empty() || cur_score < min_score) {
			min_bs = bs;
			min_score = cur_score;
		}
	}
	
	map<Node, bool>	dic;
	const vector<Node>	vs = sort_nodes();
	dic[vs[0]] = false;
	for(size_t i = 1; i < vs.size(); ++i)
		dic[vs[i]] = min_bs[i-1];
	return dic;
}

vector<InverseGraph::Edge> InverseGraph::extract_edges() const {
	vector<Edge>	edges;
	for(auto p = begin(); p != end(); ++p) {
		const Node&	v1 = p->first;
		for(auto q = p->second.begin(); q != p->second.end(); ++q) {
			const Node&	v2 = get<0>(*q);
			const int&	n1 = get<1>(*q);
			const int&	n2 = get<2>(*q);
			if(v1 < v2)
				edges.push_back(Edge(v1, v2, n1, n2));
		}
	}
	return edges;
}

int InverseGraph::calc_match_score(const vector<bool>& bs) const {
	int	num_wrong = 0;
	const auto	edges = this->extract_edges();
	const auto	dic_vs = calc_dic();
	for(auto p = edges.begin(); p != edges.end(); ++p) {
		auto	p1 = dic_vs.find(get<0>(*p));
		auto	p2 = dic_vs.find(get<1>(*p));
		const size_t	i = p1->second;
		const size_t	j = p2->second;
		const bool	b1 = i == 0 ? false : bs[i-1];
		const bool	b2 = bs[j-1];
		if(b1 == b2)
			num_wrong += get<3>(*p);
		else
			num_wrong += get<2>(*p);
	}
	return num_wrong;
}

vector<pair<double, InverseGraph::Edge>> InverseGraph::sort_edges() const {
	const vector<Edge>	edges = this->extract_edges();
	vector<pair<double, Edge>>	sorted_edges;
	for(auto p = edges.begin(); p != edges.end(); ++p) {
		const int	n1 = get<2>(*p);
		const int	n2 = get<3>(*p);
		const double	r1 = (n1 + 1) / (double)(n1 + n2 + 2);
		const double	r2 = (n2 + 1) / (double)(n1 + n2 + 2);
		const double	ratio = std::min(r1, r2);
		sorted_edges.push_back(make_pair(ratio, *p));
	}
	std::sort(sorted_edges.begin(), sorted_edges.end());
	return sorted_edges;
}

map<InverseGraph::Node, bool> InverseGraph::connect_biased_edges() const {
	const vector<pair<double, Edge>>	sorted_edges = this->sort_edges();
	// [ { node: [(index, inv?)] } ]
	vector<BoolGraph>	subgraphs;
	for(auto p = begin(); p != end(); ++p) {
		BoolGraph	graph;
		graph[p->first];
		subgraphs.push_back(graph);
	}
	
	// Add from highly biased edge
	for(auto p = sorted_edges.begin(); p != sorted_edges.end(); ++p) {
		const size_t	i = get<0>(p->second);
		const size_t	j = get<1>(p->second);
		const int		n1 = get<2>(p->second);
		const int		n2 = get<3>(p->second);
		int	k1;
		for(k1 = 0; ; ++k1) {
			const auto&	g1 = subgraphs[k1];
			if(g1.find(i) != g1.end())
				break;
		}
		int	k2;
		for(k2 = 0; ; ++k2) {
			const auto&	g2 = subgraphs[k2];
			if(g2.find(j) != g2.end())
				break;
		}
		if(k1 == k2) {
			auto&	subg1 = subgraphs[k1];
			subg1[i].push_back(make_pair(j, n1 < n2));
			subg1[j].push_back(make_pair(i, n1 < n2));
			if(!InverseGraph::is_consistent(subgraphs[k1]))
				return map<size_t, bool>();
		}
		else {
			auto&	subg1 = subgraphs[k1];
			auto&	subg2 = subgraphs[k2];
			subg1.join(subg2);
			subg1[i].push_back(make_pair(j, n1 < n2));
			subg1[j].push_back(make_pair(i, n1 < n2));
			subgraphs.erase(subgraphs.begin() + k2);
			if(subgraphs.size() == 1)
				return InverseGraph::invs(subgraphs[0]);
		}
	}
	return map<size_t, bool>();		// not coming here
}

bool InverseGraph::is_consistent(
						const map<Node, vector<pair<Node, bool>>>& graph) {
	map<Node, bool>	visited;
	Node	v0 = graph.begin()->first;
	stack<pair<Node, bool>>	stk;
	stk.push(make_pair(v0, false));
	while(!stk.empty()) {
		const Node	v = stk.top().first;
		const bool	b = stk.top().second;
		stk.pop();
		auto	iter = graph.find(v);
		for(auto p = iter->second.begin(); p != iter->second.end(); ++p) {
			const Node	v1 = p->first;
			const bool	b1_edge = p->second;
			const bool	b1 = b ^ b1_edge;
			auto	q = visited.find(v1);
			if(q != visited.end()) {
				const bool	b1_prev = q->second;
				if(b1 != b1_prev)
					return false;
			}
			else {
				visited.insert(make_pair(v1, b1));
				stk.push(make_pair(v1, b1));
			}
		}
	}
	return true;
}

map<InverseGraph::Node, bool> InverseGraph::invs(
							const map<Node, vector<pair<Node, bool>>>& graph) {
	map<Node, bool>	visited;
	const Node	v0 = graph.begin()->first;
	stack<pair<Node, bool>>	stk;
	stk.push(make_pair(v0, false));		// assume that v0 is not inverted
	while(!stk.empty()) {
		const Node	v = stk.top().first;
		const bool	b = stk.top().second;
		stk.pop();
		auto	iter = graph.find(v);
		for(auto p = iter->second.begin(); p != iter->second.end(); ++p) {
			const Node	v1 = p->first;
			const bool	b1_edge = p->second;
			const bool	b1 = b ^ b1_edge;
			if(visited.find(v1) == visited.end()) {
				visited.insert(make_pair(v1, b1));
				stk.push(make_pair(v1, b1));
			}
		}
	}
	return visited;
}

ostream& operator <<(ostream& os, const InverseGraph& graph) {
	for(auto p = graph.begin(); p != graph.end(); ++p) {
		os << p->first << " : ";
		for(auto q = p->second.begin(); q != p->second.end(); ++q)
			os << "(" << get<0>(*q) << ", " << get<1>(*q)
							<< ", " << get<2>(*q) << ") ";
		os << "\n";
	}
	return os;
}


//////////////////// BoolGraph ////////////////////

vector<BoolGraph::Node> BoolGraph::collect_nodes() const {
	vector<Node>	nodes;
	for(auto p = begin(); p != end(); ++p)
		nodes.push_back(p->first);
	return nodes;
}

vector<BoolGraph::Node> BoolGraph::neighbors(Node v0) const {
	auto	p = find(v0);
	if(p == end())
		return vector<Node>();
	else {
		vector<Node>	neighs;
		for(auto q = p->second.begin(); q != p->second.end(); ++q)
			neighs.push_back(get<0>(*q));
		return neighs;
	}
}

void BoolGraph::join(const BoolGraph& graph) {
	for(auto p = graph.begin(); p != graph.end(); ++p)
		insert(*p);
}
