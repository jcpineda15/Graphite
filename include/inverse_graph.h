#ifndef __INVERSEGRAPH
#define __INVERSEGRAPH

#include <ostream>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include "graph.h"


//////////////////// InverseGraph ////////////////////

class InverseGraph : public GraphBase,
			public 	std::map<std::size_t,
							 std::vector<std::tuple<std::size_t, int, int>>> {
public:
	typedef std::size_t	Node;
	typedef std::tuple<Node, Node, int, int>	Edge;
	
public:
	InverseGraph() { }
	~InverseGraph() { }
	
	InverseGraph *copy() const;
	
	///// virtual methods for GraphBase /////
	std::vector<Node> collect_nodes() const;
	std::vector<Node> neighbors(Node v0) const;
	
	///// non-virtual methods /////
	std::vector<bool> optimize_inversions() const;
	std::vector<InverseGraph> divide_into_connected() const;
	std::map<Node, bool> optimize_inversions_connected() const;
	
private:
	std::vector<Node> sort_nodes() const;
	std::map<Node, std::size_t> calc_dic() const;
	bool is_leaf(Node v) const;
	std::pair<InverseGraph *, std::vector<Node>> remove_acyclic_nodes() const;
	bool is_tree() const;	// assume connected
	void remove_directed_edge(Node v1, Node v2);
	void remove_edge(Node v1, Node v2);
	std::map<Node, bool> decide_inversions_tree() const;
	void add_removed_nodes(std::map<Node, bool>& dic_bs,
							const std::vector<Node>& removed_nodes) const;
	std::map<Node, bool> optimize_inversions_connected_core() const;
	std::map<Node, bool> search_all() const;
	int choice_seed() const;
	std::map<Node, bool> search_randomly() const;
	std::vector<Edge> extract_edges() const;
	int calc_match_score(const std::vector<bool>& bs) const;
	std::vector<std::pair<double, Edge>> sort_edges() const;
	std::map<Node, bool> connect_biased_edges() const;
	
public:
	static bool is_consistent(
			const std::map<Node, std::vector<std::pair<Node, bool>>>& graph);
	static std::map<Node, bool> invs(
			const std::map<Node, std::vector<std::pair<Node, bool>>>& graph);
};

std::ostream& operator <<(std::ostream& os, const InverseGraph& graph);


//////////////////// BoolGraph ////////////////////

class BoolGraph : public GraphBase,
					public std::map<std::size_t,
									std::vector<std::pair<std::size_t, bool>>> {
public:
	typedef std::size_t	Node;
	
public:
	///// virtual methods for GraphBase /////
	std::vector<Node> collect_nodes() const;
	std::vector<Node> neighbors(Node v0) const;
	
	///// non-virtual methods /////
	void join(const BoolGraph& graph);
};
#endif
