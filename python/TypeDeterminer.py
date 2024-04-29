from __future__ import annotations

# coding: utf-8
# TypeDeterminer.py

from itertools import product, count
from collections import defaultdict
from queue import PriorityQueue
from enum import Enum
from typing import Dict, Generator


#################### ParentComb ####################

class ParentComb(Enum):
	P00x00 = 0
	P00x01 = 1
	P01x01 = 2
	P00x11 = 3
	P01x11 = 4
	P11x11 = 5
	PNA = 6
	
	def is_NA(self) -> bool:
		return self == ParentComb.PNA
	
	def is_homohomo(self) -> bool:
		return self in (ParentComb.P00x00, ParentComb.P00x11, ParentComb.P11x11)
	
	def is_heterohomo(self) -> bool:
		return self in (ParentComb.P00x01, ParentComb.P01x11)
	
	def is_same_parent_genotype(self):
		return self in (ParentComb.P00x00, ParentComb.P01x01, ParentComb.P11x11)
	
	def int_gt_pair(self) -> tuple[int, int]:
		p = self.value
		if p == 0:
			return (0, 0)
		elif p < 3:
			return (p - 1, 1)
		else:
			return (p - 3, 2)


#################### TypeDeterminer ####################

class TypeDeterminer:
	def __init__(self, n: int, alpha_: float):
		self.N: int = n
		self.alpha: float = alpha_
		# どのGenotypeの組み合わせかを4進で表す
		# e.g. 0/1 x 1/1 -> 1 + 2*4 = 9
		self.memo: Dict[tuple[int, int, int], list[tuple[ParentComb, float]]] \
															= defaultdict(list)
		self.make_memo00()
		self.make_memo01()
		self.make_memo02()
		self.make_memo11()
		self.make_memo12()
		self.make_memo22()
		self.sort()
	
	@staticmethod
	def gen_errors(n, k):
		if k == 0:
			yield ()
		else:
			for n1 in range(n + 1):
				for ns in TypeDeterminer.gen_errors(n - n1, k - 1):
					yield (n1,) + ns
	
	# N/Aを含めて2割まで、N/Aを含めなければ1割まで
	@staticmethod
	def gen_error_combinations(N, k):
		for num_NA in range(N//5 + 1):
			M = min(N//10, N//5 - num_NA)
			for ns in TypeDeterminer.gen_errors(M, k - 1):
				yield (num_NA,) + ns
	
	def make_memo00(self):
		# 全部同じになるはずのパターンはその他が2割まで
		for num_NA, n1, n2 in TypeDeterminer.gen_errors(self.N//5, 3):
			n0 = self.N - n1 - n2 - num_NA
			self.memo[(n0, n1, n2)].append((ParentComb.P00x00, 0.0))
	
	def make_memo01(self):
		for num_NA, n2 in TypeDeterminer.gen_error_combinations(self.N, 2):
			M = self.N - num_NA - n2
			ps = self.binomial(M)
			for n0, p in ps:
				self.memo[(n0, M-n0, n2)].append((ParentComb.P00x01, p))
	
	def make_memo02(self):
		for num_NA, n0, n2 in TypeDeterminer.gen_errors(self.N//5, 3):
			n1 = self.N - n0 - n2 - num_NA
			self.memo[(n0, n1, n2)].append((ParentComb.P00x11, 0.0))
	
	# C++と同じ挙動を示すように大きい方が優先するように符号を反転する
	class PQ:
		def __init__(self):
			self.pq = PriorityQueue()
		
		def put(self, v):
			p, n1, n2, n3 = v
			self.pq.put((-p, -n1, -n2, -n3))
		
		def get(self):
			v = self.pq.get()
			return (-v[0], -v[1], -v[2], -v[3])
	
	# 0/1 x 0/1は難しい
	def make_memo11(self):
		# 確率が大きい状態から並べて累積が1-αを超えるまで列挙する
		for num_NA in range(self.N//5 + 1):
			M = self.N - num_NA
			pqs0 = TypeDeterminer.initialize_state(M)
			pq = TypeDeterminer.PQ()
			pq.put(pqs0)
			visited = set([pqs0])
			total_p = 0.0
			while total_p < 1.0 - self.alpha:
				p, n1, n2, n3 = pq.get()
				if n1 == n3:
					self.memo[(n1, n2, n3)].append((ParentComb.P01x01, total_p))
					total_p += p
				else:
					# n0 > n2なので、n0 < n2の分も考える
					self.memo[(n1, n2, n3)].append((ParentComb.P01x01, total_p))
					self.memo[(n3, n2, n1)].append((ParentComb.P01x01, total_p))
					total_p += p * 2
				
				neighbors = TypeDeterminer.neighbor_states((p, n1, n2, n3))
				for pqs1 in neighbors:
					if pqs1 not in visited:
						pq.put(pqs1)
						visited.add(pqs1)
	
	def make_memo12(self):
		for num_NA, n0 in TypeDeterminer.gen_error_combinations(self.N, 2):
			M = self.N - num_NA - n0
			ps = self.binomial(M)
			for n1, p in ps:
				self.memo[(n0, n1, M-n1)].append((ParentComb.P01x11, p))
	
	def make_memo22(self):
		for num_NA, n0, n1 in TypeDeterminer.gen_errors(self.N//5, 3):
			n2 = self.N - n0 - n1 - num_NA
			self.memo[(n0, n1, n2)].append((ParentComb.P11x11, 0.0))
	
	def binomial(self, M: int) -> list[tuple[int, float]]:
		p: float = 0.5**M
		ps: list[float] = [p]
		for n in range(1, M + 1):
			ps.append(ps[-1] * (M - n + 1) / n)
		
		total_p: float = 0.0
		ps2: list[tuple[int, float]] = []
		if M % 2 == 0:
			ps2.append((M//2, total_p))
			total_p += ps[M//2]
		for n1 in range((M-1)//2, -1, -1):
			n2 = M - n1
			ps2.append((n1, total_p))
			ps2.append((n2, total_p))
			total_p += ps[n1] * 2
			if total_p >= 1.0 - self.alpha:
				break
		return ps2
	
	def sort(self):
		self.memo = { ns: sorted(v, key=lambda k: k[1])
									for ns, v in self.memo.items() }
	
	def determine(self, counter: tuple[int,int,int]
									) -> list[tuple[ParentComb, float]]:
		return self.memo.get(counter, [])
	
	@staticmethod
	def genotype_probability(n0: int, n1: int, n2: int) -> float:
		# nが大きければ対数計算しなければならないかも
		n = n0 + n1 + n2
		p0 = 0.25
		p1 = 0.5
		p2 = 0.25
		p = 1.0
		for i in range(1, n0 + 1):
			p *= p0 * (n - i + 1) / i
		for i in range(1, n1 + 1):
			p *= p1 * (n - n0 - i + 1) / i
		for i in range(1, n2 + 1):
			p *= p2 * (n2 - i + 1) / i
		return p
	
	@staticmethod
	def initialize_state(M: int) -> tuple[float, int, int, int]:
		n1 = M // 2
		n0 = (M - n1 + 1) // 2	# n0 >= n2とする
		n2 = M - n0 - n1
		return TypeDeterminer.create_pqstate(n0, n1, n2)
	
	@staticmethod
	def create_pqstate(n0: int, n1: int, n2: int
								) -> tuple[float, int, int, int]:
		p = TypeDeterminer.genotype_probability(n0, n1, n2)
		return (p, n0, n1, n2)
	
	@staticmethod
	def neighbor_states(s0: tuple[float, int, int, int]
										) -> list[tuple[float, int, int, int]]:
		neighbors = []
		p, n1, n2, n3 = s0
		if n1 > 0:
			if n1-1 >= n3:
				neighbors.append(TypeDeterminer.create_pqstate(n1-1, n2+1, n3))
			if n1-1 >= n3+1:
				neighbors.append(TypeDeterminer.create_pqstate(n1-1, n2, n3+1))
		if n2 > 0:
			neighbors.append(TypeDeterminer.create_pqstate(n1+1, n2-1, n3))
			if n1 >= n3+1:
				neighbors.append(TypeDeterminer.create_pqstate(n1, n2-1, n3+1))
		if n3 > 0:
			neighbors.append(TypeDeterminer.create_pqstate(n1+1, n2, n3-1))
			neighbors.append(TypeDeterminer.create_pqstate(n1, n2+1, n3-1))
		return neighbors
