#ifndef __GRAPH
#define __GRAPH

#include <vector>
#include <map>


//////////////////// GraphBase ////////////////////

class GraphBase {
	using Index = std::size_t;
	using Node = std::size_t;
	using Edge = std::pair<Node, Node>;
	
	struct EdgeGenerator {
		const GraphBase& graph;
		const std::vector<Node>	nodes;
		std::vector<Node>::const_iterator	p;
		std::vector<Node>	neighs;
		std::vector<Node>::const_iterator	q;
		
		EdgeGenerator(const GraphBase& g);
		
		bool ends() const;
		bool is_effective() const;
		void proceed();
		void find_next_edge();
		Edge next();
	};
	
public:
	GraphBase() { }
	virtual ~GraphBase() { }
	
	std::vector<std::vector<Node>> divide_nodes_into_connected() const;
	
	virtual std::vector<Node> collect_nodes() const = 0;
	virtual std::vector<Node> neighbors(Node v0) const = 0;
	
public:
	static std::size_t count_groups(const std::map<Node, std::size_t>& groups);
};


//////////////////// WeightedGraphBase ////////////////////

template<typename T>
class WeightedGraphBase : public GraphBase {
public:
	using Node = std::size_t;
	using Edge = std::tuple<Node, Node, T>;
	
	struct LessWeight {
		bool operator ()(const Edge& e1, const Edge& e2) const {
			return std::get<2>(e1) < std::get<2>(e2);
		}
	};
	
public:
	virtual std::vector<Edge> collect_weighted_edges() const = 0;
};

#endif
