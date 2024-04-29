from __future__ import annotations

# coding: utf-8
# VCFHeteroHomo.py

from collections import defaultdict, Counter
from math import log10
from typing import List, Tuple, Optional, IO, Dict, Iterator
import random
import time

from VCFFamily import *
from VCFImpFamily import FillType, VCFImpFamilyRecord
from TypeDeterminer import ParentComb
from Map import *
import Imputer
from option import *
from inverse_graph import *
from graph import Node
from invgraph import InvGraph
from common import divide_graph_into_connected


#################### VCFHeteroHomoRecord ####################

class VCFHeteroHomoRecord(VCFImpFamilyRecord):
	def __init__(self, v: list[str], samples: list[str],
							i: int, parents_wrong_type: str, pair: ParentComb):
		super().__init__(v, samples, i, parents_wrong_type, pair)
		self.which_comes_from: list[int] = [-1] * self.num_progenies()
	
	def is_mat_hetero(self) -> bool:
		return self.mat_int_gt() == 1
	
	def is_pat_hetero(self) -> bool:
		return self.pat_int_gt() == 1
	
	def is_imputable(self) -> bool:
		return self.parents_wrong_type == 'Right'
	
	def get_fill_type(self) -> FillType:
		if self.is_imputable():
			return FillType.MAT if self.is_mat_hetero() else FillType.PAT
		else:
			return FillType.IMPUTABLE
	
	def genotypes_from_hetero_parent(self) -> list[int]:
		def encode(gt, homo_parent):
			gt_ = gt if homo_parent == 0 else gt - 1
			if gt_ not in (0, 1):
				return -1
			else:
				return gt_
		
		gts = self.get_int_gts()
		homo = self.pat_int_gt() if self.is_mat_hetero() else self.mat_int_gt()
		return [ encode(gt, homo) for gt in gts[2:] ]
	
	def set_haplo(self, h: int):
		hetero_col = 9 if self.mat_int_gt() == 1 else 10
		homo_col = 10 if hetero_col == 9 else 9
		self.v[hetero_col] = '0|1' if h == 0 else '1|0'
		self.v[homo_col] = '0|0' if self.v[homo_col][0] == '0' else '1|1'
		
		gts = self.genotypes_from_hetero_parent()
		for i, gt in enumerate(gts):
			self.which_comes_from[i] = -1 if gt == -1 else 0 if gt == h else 1
	
	def set_int_gt_by_which_comes_from(self, ws: list[int]):
		self.which_comes_from = ws
		mat_gt: str = self.v[9]
		pat_gt: str = self.v[10]
		if self.is_mat_hetero():
			pat_int_gt = int(pat_gt[0])
			for i in range(self.num_progenies()):
				self.v[i+11] = '%s|%d' % (mat_gt[ws[i]*2], pat_int_gt)
		else:
			mat_int_gt = int(mat_gt[0])
			for i in range(self.num_progenies()):
				self.v[i+11] = '%d|%s' % (mat_int_gt, pat_gt[ws[i]*2])


#################### VCFHeteroHomo ####################

class VCFHeteroHomo(VCFBase, VCFSmallBase, VCFFamilyBase, VCFMeasurable):
	def __init__(self, header: list[list[str]],
					records: list[VCFHeteroHomoRecord], map_: Map):
		self.records = records
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		VCFFamilyBase.__init__(self)
		VCFMeasurable.__init__(self, map_)
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def is_mat_hetero(self):
		return self.records[0].is_mat_hetero()
	
	def make_graph(self, max_dist: float) -> InvGraph:
		# ベータ分布を使うので、distはfloat
		def distance(gts1: list[int], gts2: list[int]) -> tuple[float, bool]:
			def dist_with_NA(right: int, counter_NA: int) -> float:
				diff = N - right - counter_NA
				diff_ratio = (diff + 1.0) / (right + diff + 2.0)
				return diff + counter_NA * diff_ratio
			
			N = len(gts1)
			counter_right = 0
			counter_diff = 0
			counter_NA = 0
			for int_gt1, int_gt2 in zip(gts1, gts2):
				if int_gt1 == -1 or int_gt2 == -1:
					counter_NA += 1
				elif int_gt1 == int_gt2:
					counter_right += 1
				else:
					counter_diff += 1
			
			if counter_right >= counter_diff:
				return (dist_with_NA(counter_right, counter_NA), False)
			else:
				return (dist_with_NA(counter_diff, counter_NA), True)
		
		L = len(self)
		graph = InvGraph()
		for k in range(L):
			graph[Node(k)] = []
		gtss = [ r.genotypes_from_hetero_parent() for r in self.records ]
		cMs = [ self.cM(record.pos()) for record in self.records ]
		# C++版ではkで並列化したい
		for k in range(L):
			for l in range(k+1, L):
				cM = cMs[l] - cMs[k]
				if cM > 10.0:	# 10cM以上離れていたら繋がりを見ない
					break
				d, b = distance(gtss[k], gtss[l])
				if d <= max_dist:
					graph[k].append((l, d, b))
					graph[l].append((k, d, b))
		return graph
	
	# haplotype1の各レコードのGenotypeが0なのか1なのか
	def make_parent_haplotypes(self, graph: InvGraph) -> list[int]:
		def is_visited(v: int) -> bool:
			return haplo[v] != -1
		
		tree: InvGraph = graph.minimum_spanning_tree()
		vs = sorted(tree.keys())
		haplo = self.create_haplotype(vs[0], tree)
		return [ h for h in haplo if h != -1 ]
	
	def create_haplotype(self, v0: Node, tree: InvGraph) -> list[int]:
		L = len(self)
		haplo: list[int] = [-1] * L
		haplo[v0] = 0
		edges = tree.walk(v0)
		for v1, v2, inv in edges:
			if haplo[v2] == -1:
				if inv:
					haplo[v2] = 1 if haplo[v1] == 0 else 0
				else:
					haplo[v2] = haplo[v1]
			elif haplo[v1] == -1:
				if inv:
					haplo[v1] = 1 if haplo[v2] == 0 else 0
				else:
					haplo[v1] = haplo[v2]
		return haplo
	
	def make_subvcf(self, graph: InvGraph) -> VCFHeteroHomo:
		haplo = self.make_parent_haplotypes(graph)
		indices = sorted(graph.keys())
		records = [ self.records[i] for i in indices ]
		for record, h in zip(records, haplo):
			record.set_haplo(h)
		
		return VCFHeteroHomo(self.header, records, self.map)
	
	def determine_haplotype(self, option: OptionImpute
					) -> tuple[list[VCFHeteroHomo], list[VCFHeteroHomoRecord]]:
		max_dist = min(option.max_dist, (self.num_samples() - 2) * 0.1)
		graph = self.make_graph(max_dist)
		subgraphs = graph.divide_into_connected()
		
		# 小さなgraphのmarkerは省く
		# 小さなgraphしかなければ、その中でも一番大きなgraphにする
		gs = [ g for g in subgraphs if len(g) >= option.min_graph ]
		if not gs:
			max_index, g = max(enumerate(subgraphs), key=lambda v: len(v[1]))
			gs = [g]
			# gに使われたrecordは入れない
			unused_records = [ self.records[j]
								for i, g in enumerate(subgraphs)
								if i != max_index and len(g) < option.min_graph
								for j in g.keys() ]
		else:
			unused_records = [ self.records[i]
							for g in subgraphs if len(g) < option.min_graph
							for i in g.keys() ]
		
		subvcfs = list(map(self.make_subvcf, gs))
		return (subvcfs, unused_records)
	
	def make_seq(self, i):
		def convert(gt):
			if gt == 0 or gt == 1:
				return str(gt)
			else:
				return 'N'
		
		return ''.join(convert(record.which_comes_from[i])
									for record in self.records)
	
	def impute_each_sample_seq(self, i: int, cMs: list[float], min_c: float):
		def is_all_same_without_N(seq):
			prev_c = '-'
			for c in seq:
				if c == 'N':
					continue
				elif prev_c == '-':
					prev_c = c
				elif c != prev_c:
					return False
			else:
				return True
		
		def create_same_color_string(seq):
			if all(c == 'N' for c in seq):
				c_base = '0'
			else:
				c_base = next(c for c in seq if c != 'N')
			return c_base * len(seq)
		
		seq = self.make_seq(i)
		if is_all_same_without_N(seq):
			return create_same_color_string(seq)
		
		hidden_states = ['0', '1']
		states = ['0', '1', 'N']
		hidden_seq = Imputer.impute(seq, hidden_states, states, cMs)
		painted_seq = Imputer.paint(hidden_seq, cMs, min_c)
		return painted_seq
	
	def impute_each(self, option: OptionImpute):
		cMs = [ self.cM(record.pos()) for record in self.records ]
		imputed_seqs = [
				self.impute_each_sample_seq(i, cMs, option.min_crossover)
								for i in range(self.num_samples() - 2) ]
		for k, record in enumerate(self.records):
			ws = [ int(seq[k]) for seq in imputed_seqs ]
			record.set_int_gt_by_which_comes_from(ws)
	
	def create_option(self) -> OptionImpute:
		num = max(2, len(self))
		total_cM = self.cM(self.records[-1].pos())
		max_dist = max(4, int(total_cM * self.num_progenies()
											/ num * log10(num) * 2.5 * 0.01))
		return OptionImpute(max_dist, 20, 5, 1.0)
	
	def impute(self) -> tuple[list[VCFHeteroHomo], list[VCFHeteroHomoRecord]]:
		if not self.records:
			return ([self], [])
		option: OptionImpute = self.create_option()
		vcfs, unused_records = self.determine_haplotype(option)
		for vcf in vcfs:
			vcf.impute_each(option)
		for record in unused_records:
			record.enable_modification()
		return (vcfs, unused_records)
	
	# 共通のヘテロ親はどれだけマッチしているか
	def match(self, other: VCFHeteroHomo) -> tuple[int, int]:
		if len(self) == 0 or len(other) == 0:
			return (0, 0)
		
		hetero_col1 = 9 if self.is_mat_hetero() else 10
		hetero_col2 = 9 if other.is_mat_hetero() else 10
		L, M = len(self), len(other)
		num_match, num_unmatch = 0, 0
		k, l = 0, 0
		while k < L and l < M:
			record1 = self.records[k]
			record2 = other.records[l]
			if record1.pos() == record2.pos():
				if record1.v[hetero_col1] == record2.v[hetero_col2]:
					num_match += 1
				else:
					num_unmatch += 1
			if record1.pos() <= record2.pos():
				k += 1
			if record1.pos() >= record2.pos():
				l += 1
		
		return (num_match, num_unmatch)
	
	def inverse_hetero_parent_phases(self):
		hetero_col = 9 if self.is_mat_hetero() else 10
		for record in self.records:
			if record.v[hetero_col] == '0|1':
				record.v[hetero_col] = '1|0'
			else:
				record.v[hetero_col] = '0|1'
	
	# Create a pair of VCFHeteroHomo for each Family and store it for each parent
	@staticmethod
	def make_VCFHeteroHomo(records: list[VCFHeteroHomoRecord],
							header: list[list[str]], geno_map: Map
									) -> tuple[VCFHeteroHomo, VCFHeteroHomo,
													list[VCFHeteroHomoRecord]]:
		unused_records: list[VCFHeteroHomoRecord] = []
		heho_mat_records = [ r for r in records
								if r.is_imputable() and r.is_mat_hetero() ]
		heho_pat_records = [ r for r in records
								if r.is_imputable() and not r.is_mat_hetero() ]
		# TODO: classifyのときにunusedに入れておけばいいのでは？
		unused_records = [ r for r in records if not r.is_imputable() ]
		vcf_mat = VCFHeteroHomo(header, heho_mat_records, geno_map)
		vcf_pat = VCFHeteroHomo(header, heho_pat_records, geno_map)
		return (vcf_mat, vcf_pat, unused_records)
	
	@staticmethod
	def impute_vcfs(records: list[VCFHeteroHomoRecord],
					header: list[list[str]], geno_map: Map
					) -> tuple[list[VCFHeteroHomo], list[VCFHeteroHomoRecord]]:
		vcf_mat, vcf_pat, unused = VCFHeteroHomo.make_VCFHeteroHomo(
													records, header, geno_map)
		vcfs_mat, unused_records1 = vcf_mat.impute()
		vcfs_pat, unused_records2 = vcf_pat.impute()
		vcfs = []
		vcfs.extend(vcfs_mat)
		vcfs.extend(vcfs_pat)
		unused.extend(unused_records1)
		unused.extend(unused_records2)
		return (vcfs, unused)
	
	@staticmethod
	def inverse_phases(vcfs: list[VCFHeteroHomo]):
		graph = VCFHeteroHomo.make_vcf_graph(vcfs)
		if not graph:
			return
		
		bs: tuple[bool, ...] = graph.optimize_inversions()
		for vcf, b in zip(vcfs, bs):
			if b:
				vcf.inverse_hetero_parent_phases()
	
	@staticmethod
	def make_vcf_graph(vcfs: list[VCFHeteroHomo]) -> InverseGraph:
		N = len(vcfs)
		graph: InverseGraph = InverseGraph()
		for v in range(N):
			graph[Node(v)] = []
		for i, j in combinations(range(N), 2):
			num_match, num_unmatch = vcfs[i].match(vcfs[j])
			if num_match != 0 or num_unmatch != 0:
				graph[i].append((j, num_match, num_unmatch))
				graph[j].append((i, num_match, num_unmatch))
		return graph

__all__ = ['VCFHeteroHomoRecord', 'VCFHeteroHomo']
