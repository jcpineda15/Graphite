#ifndef __INVGRAPH
#define __INVGRAPH

#include <vector>
#include <map>
#include "graph.h"


//////////////////// InvGraph ////////////////////

class InvGraph : public WeightedGraphBase<double>,
			public std::map<std::size_t,
						std::vector<std::tuple<std::size_t, double, bool>>> {
public:
	using Node = std::size_t;
	using Edge = std::tuple<Node, Node, bool>;
	
public:
	InvGraph() { }
	~InvGraph() { }
	
	///// virtual methods for GraphBase /////
	std::vector<Node> collect_nodes() const;
	std::vector<Node> neighbors(Node v0) const;
	
	///// virtual methods for WeightedGraphBase /////
	std::vector<WeightedGraphBase<double>::Edge> collect_weighted_edges() const;
	
	///// non-virtual methods /////
	InvGraph minimum_spanning_tree() const;
	std::vector<InvGraph> divide_into_connected() const;
	std::vector<Edge> walk(Node v0) const;
	
public:
	static std::size_t count_groups(const std::map<Node, std::size_t>& groups);
};
#endif
