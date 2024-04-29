from __future__ import annotations

# coding: utf-8
# UnionFind.py

from typing import NewType, Tuple, Iterator


#################### Type ####################

Node = NewType('Node', int)


#################### UnionFind ####################

class UnionFind:
	def __init__(self, nodes: list[Node]):
		self.parents: dict[Node, Node] = { v: v for v in nodes }
		self.heights: dict[Node, int] = { v: 1 for v in nodes }
	
	def join(self, v1: Node, v2: Node):
		r1 = self.root(v1)
		r2 = self.root(v2)
		h1 = self.heights[r1]
		h2 = self.heights[r2]
		if h1 <= h2:
			self.parents[r1] = r2
			self.heights[r2] = max(self.heights[r2], self.heights[r1]+1)
		else:
			self.parents[r2] = r1
			self.heights[r1] = max(self.heights[r1], self.heights[r2]+1)
	
	def root(self, v0: Node) -> Node:
		v = v0
		while True:
			p = self.parents[v]
			if p == v:
				return p
			v = p
	
	def divide_into_trees(self) -> dict[Node, int]:
		groups: dict[Node, int] = { }
		group_index = 0
		for v0, p in self.parents.items():
			if v0 in groups:
				continue
			
			# find root or defined node
			path: list[Node] = [v0]
			while True:
				v = path[-1]
				if v in groups:
					break
				
				p1 = self.parents[v]
				if p1 == v:		# root
					groups[v] = group_index
					group_index += 1
					break
				else:
					path.append(p1)
			
			for v in path[:-1]:
				groups[v] = groups[path[-1]]
		
		return groups


#################### UnionFindSerial ####################

class UnionFindSerial:
	def __init__(self, N: int):
		self.parents: list[int] = list(range(N))
		self.heights: list[int] = [1] * (N)
	
	def join(self, v1: int, v2: int):
		r1 = self.root(v1)
		r2 = self.root(v2)
		h1 = self.heights[r1]
		h2 = self.heights[r2]
		if h1 <= h2:
			self.parents[r1] = r2
			self.heights[r2] = max(self.heights[r2], self.heights[r1]+1)
		else:
			self.parents[r2] = r1
			self.heights[r1] = max(self.heights[r1], self.heights[r2]+1)
	
	def root(self, v: int) -> int:
		while self.parents[v] != v:
			v = self.parents[v]
		return v
	
	def __str__(self) -> str:
		def f(v: int) -> str:
			s = str(v)
			for child in children[v]:
				for line in f(child).split('\n'):
					s += '\n  ' + line
			return s
		
		N = len(self.parents)
		children: list[list[int]] = [ [] for v in range(N) ]
		for v in range(N):
			if self.parents[v] != v:
				children[self.parents[v]].append(v)
		
		return '\n'.join(f(v) for v in range(N) if self.parents[v] == v)

__all__ = ['UnionFind', 'UnionFindSerial', 'Node']
