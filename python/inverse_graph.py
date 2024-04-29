from __future__ import annotations

# coding: utf-8
# inverse_graph.py
# Nodeを反転すべきかどうかを決めるグラフ

from itertools import count, product
from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator
import random
import time

from graph import GraphBase, Node


#################### Type ####################

Edge = Tuple[Node, Node, int, int]
Value = Tuple[Node, int, int]


def is_consistent(graph: Dict[int, list[tuple[int, bool]]]) -> bool:
	visited: Dict[int, bool] = { }
	v0 = next(v for v in graph.keys())
	stk = [(v0, False)]		# v0は反転しないとする
	while stk:
		v, b = stk.pop()
		for v1, b1_edge in graph[v]:
			b1 = b ^ b1_edge
			if v1 in visited:
				b1_prev = visited[v1]
				if b1 != b1_prev:
					return False
			else:
				visited[v1] = b1
				stk.append((v1, b1))
	return True
	
def invs(graph: Dict[Node, list[tuple[Node, bool]]]) -> Dict[Node, bool]:
	visited: Dict[Node, bool] = { }
	v0 = min(v for v in graph.keys())
	stk = [(v0, False)]		# v0は反転しないとする
	while stk:
		v, b = stk.pop()
		for v1, b1_edge in graph[v]:
			b1 = b ^ b1_edge
			if v1 not in visited:
				visited[v1] = b1
				stk.append((v1, b1))
	return visited


#################### InverseGraph ####################

class InverseGraph(dict, GraphBase):
	def __init__(self, *args, **kwargs):
		dict.__init__(self, *args, **kwargs)
		GraphBase.__init__(self)
	
	def __setitem__(self, key: Node, value: list[Value]):
		super().__setitem__(key, value)
	
	def copy(self) -> InverseGraph:
		graph = InverseGraph()
		for v, vs in self.items():
			graph[v] = vs[:]
		return graph
	
	##### virtual methods for GraphBase #####
	def generate_nodes(self) -> Iterator[Node]:
		for v in self.keys():
			yield v
	
	def neighbors(self, v0: Node) -> list[Node]:
		value = self.get(v0, [])
		return [ v for v, n1, n2 in value ]
	
	##### non-virtual methods #####
	def optimize_inversions(self) -> tuple[bool, ...]:
		subgraphs = self.divide_graph_into_connected()
		bs = [False] * len(self)
		for subg in subgraphs:
			dic_bs = subg.optimize_inversions_connected()
			for i, b in dic_bs.items():
				bs[i] = b
		return tuple(bs)
	
	def divide_graph_into_connected(self) -> list[InverseGraph]:
		vss = self.divide_nodes_into_connected()
		subgraphs: list[InverseGraph] = []
		for vs in vss:
			subg = InverseGraph()
			for v in vs:
				subg[v] = self.get(v, [])
			subgraphs.append(subg)
		return subgraphs
	
	def is_leaf(self, v: Node) -> bool:
		return v in self and len(self[v]) == 1
	
	def remove_directed_edge(self, v1: Node, v2: Node):
		neighs1 = self.neighbors(v1)
		if len(neighs1) == 1:
			del self[v1]
		else:
			vs2 = self[v1]
			for i, (v, n1, n2) in enumerate(vs2):
				if v == v2:
					vs2.pop(i)
					break
			else:
				raise Exception('error in InverseGraph.remove_directed_edge')
	
	def remove_edge(self, v1: Node, v2: Node):
		self.remove_directed_edge(v1, v2)
		self.remove_directed_edge(v2, v1)
	
	def remove_acyclic_nodes(self) -> tuple[InverseGraph, list[Node]]:
		graph = self.copy()
		removed_nodes: list[Node] = []
		while True:
			prev_num_removed_nodes = len(removed_nodes)
			vs = list(graph.keys())
			for v in vs:
				if graph.is_leaf(v):
					neighs = graph.neighbors(v)
					v1 = neighs[0]
					graph.remove_edge(v, v1)
					removed_nodes.append(v)
			
			if len(removed_nodes) == prev_num_removed_nodes:
				break
		return (graph, removed_nodes)
	
	# assume connected
	def is_tree(self) -> bool:
		num_edges = sum(len(vs) for vs in self.values())
		return num_edges == (len(self) - 1) * 2
	
	def decide_inversions_tree(self) -> Dict[Node, bool]:
		dic_bs = { }
		v0 = next(v for v in sorted(self.keys()))
		visited = set([v0])
		stk = [(v0, False)]
		dic_bs[v0] = False
		while stk:
			(v, b) = stk.pop()
			for v1, n1, n2 in self[v]:
				b1 = b ^ (n1 < n2)
				if v1 not in visited:
					dic_bs[v1] = b1
					stk.append((v1, b1))
					visited.add(v1)
		return dic_bs
	
	# Assume that self is a connected component
	def optimize_inversions_connected_core(self) -> Dict[Node, bool]:
		if len(self) <= 20:
			return self.search_all()
		
		# Finalize the edge from the one with the largest true/false bias
		# If there are no inconsistencies by the time the graphs are all connected, it's OK.
		dic_bs = self.connect_biased_edges()
		if dic_bs:
			return dic_bs
		
		# If that doesn't work, use random.
		# I really want to use SA.
		return self.search_randomly()
	
	def add_removed_nodes(self, dic_bs: Dict[Node, bool],
								removed_nodes: list[Node]):
		for v1 in reversed(removed_nodes):
			for v2, n1, n2 in self[v1]:
				if v2 in dic_bs:
					b2 = dic_bs[v2]
					b1 = b2 ^ (n1 < n2)
					dic_bs[v1] = b1
					break
	
	# Assume that self is a connected component
	def optimize_inversions_connected(self) -> Dict[Node, bool]:
		if self.is_tree():
			return self.decide_inversions_tree()
		
		graph, removed_nodes = self.remove_acyclic_nodes()
		dic_bs = graph.optimize_inversions_connected_core()
		self.add_removed_nodes(dic_bs, removed_nodes)
		return dic_bs
	
	# Brute force
	def search_all(self) -> Dict[Node, bool]:
		if len(self) == 1:
			return { list(self.keys())[0]: False }
		
		bs = min(product((False, True), repeat=len(self)-1),
										key=lambda bs: self.score(bs))
		vs = sorted(self.keys())
		return dict(zip(vs, (False,) + bs))
	
	def search_randomly(self) -> Dict[Node, bool]:
		bs = min((tuple(random.choices((False, True), k=len(self)-1))
												for _ in range(2**19)),
												key=lambda bs: self.score(bs))
		vs = sorted(self.keys())
		return dict(zip(vs, (False,) + bs))
	
	def generate_biased_edges(self) -> Iterator[Edge]:
		for i, v in self.items():
			for j, n1, n2 in v:
				if i > j:
					continue
				yield (i, j, n1, n2)
	
	def score(self, bs: tuple[bool, ...]) -> int:
		dic: Dict[int, int] = dict(zip(sorted(self.keys()), count()))
		num_wrong: int = 0
		for v1, v2, n1, n2 in self.generate_biased_edges():
			i = dic[v1]
			j = dic[v2]
			b1 = False if i == 0 else bs[i-1]
			b2 = bs[j-1]
			if b1 == b2:
				num_wrong += n2
			else:
				num_wrong += n1
		return num_wrong
	
	def connect_biased_edges(self) -> Dict[Node, bool]:
		def ratio(v: Edge) -> float:
			i, j, n1, n2 = v
			return min((n1 + 1) / (n1 + n2 + 2), (n2 + 1) / (n1 + n2 + 2))
		
		def join(graph1: Dict[int, list[tuple[int, bool]]],
							graph2: Dict[int, list[tuple[int, bool]]]):
			for v1, vs in graph2.items():
				graph1[v1] = vs
		
		# これってminimum spanning treeなのでは
		sorted_edges = sorted(((ratio(e), e)
								for e in self.generate_biased_edges()))
		subgraphs: list[BoolGraph] = []
		for v in self.keys():
			graph = BoolGraph()
			graph[v] = []
			subgraphs.append(graph)
		
		for r, (i, j, n1, n2) in sorted_edges:
			k1 = next(k for k, g in enumerate(subgraphs) if i in g.keys())
			k2 = next(k for k, g in enumerate(subgraphs) if j in g.keys())
			if k1 == k2:
				g1 = subgraphs[k1]
				g1[i].append((j, n1 < n2))
				g1[j].append((i, n1 < n2))
				if not is_consistent(g1):
					return { }
			else:
				subg1 = subgraphs[k1]
				subg2 = subgraphs.pop(k2)
				join(subg1, subg2)
				subg1[i].append((j, n1 < n2))
				subg1[j].append((i, n1 < n2))
				if len(subgraphs) == 1:
					return invs(subgraphs[0])
		else:
			return { }	# not come here
	
	@staticmethod
	def convert(graph: list[list[tuple[int, int, int]]]) -> InverseGraph:
		return InverseGraph(dict(enumerate(graph)))


#################### BoolGraph ####################

class BoolGraph(dict, GraphBase):
	def __init__(self, *args, **kwargs):
		dict.__init__(self, *args, **kwargs)
		GraphBase.__init__(self)
	
	def __setitem__(self, key: Node, value: list[Value]):
		super().__setitem__(key, value)
	
	##### virtual methods for GraphBase #####
	def generate_nodes(self) -> Iterator[Node]:
		for v in self.keys():
			yield v
	
	def neighbors(self, v0: Node) -> list[Node]:
		value = self.get(v0, [])
		return [ v for v, b in value ]
	
	##### non-virtual methods #####
	def join(self, graph: BoolGraph):
		for u, vs in graph.items():
			self.__setitem__(u, vs)


__all__ = ['InverseGraph']
