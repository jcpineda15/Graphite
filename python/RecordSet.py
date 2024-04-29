# coding: utf-8
# RecordSet.py
# 自分自身と前後のヘテロ×ホモ親のセット

from __future__ import annotations
from typing import Optional, Iterator
from itertools import product
from math import log

from VCF import VCFRecord
from VCFFillableRecord import VCFFillableRecord
from TypeDeterminer import ParentComb
from Genotype import Genotype
from common import unique_list


#################### RecordSet ####################

class RecordSet:
	def __init__(self, r: Optional[VCFFillableRecord],
			pm: Optional[VCFFillableRecord], nm: Optional[VCFFillableRecord],
			pp: Optional[VCFFillableRecord], np: Optional[VCFFillableRecord]):
		self.record: Optional[VCFFillableRecord] = r
		self.prev_mat_record: Optional[VCFFillableRecord] = pm
		self.next_mat_record: Optional[VCFFillableRecord] = nm
		self.prev_pat_record: Optional[VCFFillableRecord] = pp
		self.next_pat_record: Optional[VCFFillableRecord] = np
	
	def records(self) -> list[Optional[VCFFillableRecord]]:
		return [self.record, self.prev_mat_record, self.next_mat_record,
							 self.prev_pat_record, self.next_pat_record]
	
	def prev_mat_from(self, i: int) -> int:
		if self.prev_mat_record is None:
			return 0
		return self.prev_mat_record.from_which_chrom(i, True)
	
	def next_mat_from(self, i: int) -> int:
		if self.next_mat_record is None:
			return 0
		return self.next_mat_record.from_which_chrom(i, True)
	
	def prev_pat_from(self, i: int) -> int:
		if self.prev_pat_record is None:
			return 0
		return self.prev_pat_record.from_which_chrom(i, False)
	
	def next_pat_from(self, i: int) -> int:
		if self.next_pat_record is None:
			return 0
		return self.next_pat_record.from_which_chrom(i, False)
	
	def gt_each(self, i: int, record: Optional[VCFFillableRecord]) -> str:
		return './.' if record is None else record.v[i+9]
	
	def gt(self, i: int) -> str:
		return self.gt_each(i, self.record)
	
	def gts(self, i: int) -> tuple[str, str, str, str, str]:
		v = [ self.gt_each(i, r) for r in self.records() ]
		return (v[0], v[1], v[2], v[3], v[4])
	
	def is_prev_nearer(self, is_mat: bool) -> bool:
		if self.record is None:
			return False
		elif is_mat:
			# どちらかがNoneの場合、ここには来ないので、適当に処理する
			if self.prev_mat_record is None or self.next_mat_record is None:
				return False
			else:
				return (self.prev_mat_record.pos() + self.next_mat_record.pos()
														> self.record.pos()*2)
		else:
			if self.prev_pat_record is None or self.next_pat_record is None:
				return False
			else:
				return (self.prev_pat_record.pos() + self.next_pat_record.pos()
														> self.record.pos()*2)
	
	def __select_phasing(self, candidates: list[tuple[int, int]]
												) -> tuple[int, int]:
		if len(candidates) == 1 or self.record is None:
			return candidates[0]
		
		mat_int_gt = self.record.mat_int_gt()
		pat_int_gt = self.record.pat_int_gt()
		
		v: list[tuple[float, int, int]] = []
		for mat_phasing, pat_phasing in candidates:
			mat_int_gt1 = (mat_phasing >> 1) + (mat_phasing & 1)
			pat_int_gt1 = (pat_phasing >> 1) + (pat_phasing & 1)
			score = (abs(mat_int_gt1 - mat_int_gt) +
					 abs(pat_int_gt1 - pat_int_gt))
			v.append((score, mat_phasing, pat_phasing))
		
		# これだと、同じスコアならphasingが小さい方から選んでいる
		# 本当はランダム的に選びたい
		_, mat_phasing, pat_phasing = min(v)
		return (mat_phasing, pat_phasing)
	
	def determine_phasing_core(self, lls: list[tuple[float, int, int]]
														) -> tuple[int, int]:
		# 最大に近いllを集める
		candidates = []
		max_ll = lls[-1][0]
		for ll, mat_phasing, pat_phasing in reversed(lls):
			if ll > max_ll - 1e-9:
				candidates.append((mat_phasing, pat_phasing))
			else:
				break
		return self.__select_phasing(candidates)
	
	def possible_phasings(self) -> list[tuple[int, int]]:
		# pair 0: 0x0, 1: 0x1, 2: 0x2, 3: 1x1, 4: 1x2, 5: 2x2
		if self.record is None:
			return []
		elif self.record.pair == ParentComb.P00x00:
			return [(0, 0)]
		elif self.record.pair == ParentComb.P00x01:
			return [(0, 1), (0, 2), (1, 0), (2, 0),
					(1, 1), (1, 2), (2, 1), (2, 2)]
		elif self.record.pair == ParentComb.P01x01:
			return [(1, 1), (1, 2), (2, 1), (2, 2),
					(0, 1), (0, 2), (1, 0), (2, 0),
					(1, 3), (2, 3), (3, 1), (3, 2)]
		elif self.record.pair == ParentComb.P00x11:
			return [(0, 3), (3, 0)]
		elif self.record.pair == ParentComb.P01x11:
			return [(1, 3), (2, 3), (3, 1), (3, 2),
					(1, 1), (1, 2), (2, 1), (2, 2)]
		elif self.record.pair == ParentComb.P11x11:
			return [(3, 3)]
		else:
			# all
			return list(product(range(4), range(4)))
	
	# phasingされている前提
	def from_which_chrom(self, gt: str, record: Optional[VCFRecord],
														mat: bool) -> int:
		if record is None:
			return 0
		
		i = 0 if mat else 1
		parent_gt = record.v[9+i]
		return 1 if parent_gt[0] == gt[i*2] else 2
	
	def from_which_chrom_prev_mat(self, gt) -> int:
		return self.from_which_chrom(gt, self.prev_mat_record, True)
	
	def from_which_chrom_next_mat(self, gt) -> int:
		return self.from_which_chrom(gt, self.next_mat_record, True)
	
	def from_which_chrom_prev_pat(self, gt) -> int:
		return self.from_which_chrom(gt, self.prev_pat_record, False)
	
	def from_which_chrom_next_pat(self, gt) -> int:
		return self.from_which_chrom(gt, self.next_pat_record, False)
	
	def gen_gts(self) -> Iterator[tuple[str, str, str, str, str]]:
		if self.record is None:
			return
		
		for c in range(11, len(self.record.v)):
			gts = [ r.v[c] if r else '' for r in self.records() ]
			yield (gts[0], gts[1], gts[2], gts[3], gts[4])
	
	def is_mat_prev_near(self) -> bool:
		if (self.record is None or self.prev_mat_record is None or
										self.next_mat_record is None):
			return False
		return (self.record.pos() * 2 <
					self.prev_mat_record.pos() + self.next_mat_record.pos())
	
	def is_pat_prev_near(self) -> bool:
		if (self.record is None or self.prev_pat_record is None or
										self.next_pat_record is None):
			return False
		return (self.record.pos() * 2 <
					self.prev_pat_record.pos() + self.next_pat_record.pos())
	
	def near_mat_from(self, i: int) -> int:
		return (self.prev_mat_from(i) if self.is_mat_prev_near()
										else self.next_mat_from(i))
	
	def near_pat_from(self, i: int) -> int:
		return (self.prev_pat_from(i) if self.is_pat_prev_near()
										else self.next_pat_from(i))
	
	# [(mat_from, pat_from)] -> (mat_from, pat_from)
	def select_nearest_froms(self, pairs: list[tuple[int, int]],
											i: int) -> tuple[int, int]:
		if len(pairs) == 4:		# N/Aのときにまれにあり得る
			return (self.near_mat_from(i), self.near_pat_from(i))
		elif pairs[0][0] == pairs[1][0]:		# matが同じ
			if (self.record is None or self.prev_pat_record is None or
											self.next_pat_record is None):
				return (0, 0)
			elif self.is_pat_prev_near():
				return (pairs[0][0], self.prev_pat_from(i))
			else:
				return (pairs[0][0], self.next_pat_from(i))
		elif pairs[0][1] == pairs[1][1]:	# patが同じ
			if (self.record is None or self.prev_mat_record is None or
											self.next_mat_record is None):
				return (0, 0)
			elif self.is_mat_prev_near():
				return (self.prev_mat_from(i), pairs[0][1])
			else:
				return (self.next_mat_from(i), pairs[0][1])
		else:	# 両親とも乗り換えている（滅多にない）
			return (self.near_mat_from(i), self.near_pat_from(i))
	
	def select_pair(self, pairs: list[tuple[int, int]], i: int,
							selected: bool = False) -> tuple[int, int]:
		def sum_gt(gt: str) -> int:
			try:
				return int(gt[0]) + int(gt[2])
			except ValueError:
				return -1
		
		record = self.record
		if record is None:	# for mypy
			return (0, 0)
		gt = record.get_gt(i)
		if not pairs:
			return (0, 0)
		elif len(pairs) == 1:
			return pairs[0]
		elif not Genotype.is_valid(gt, record.mat_int_gt(),
											record.pat_int_gt()):
			return self.select_nearest_froms(pairs, i)
		elif selected:
			return self.select_nearest_froms(pairs, i)
		else:
			new_pairs = [ v for v in pairs
						if sum_gt(gt) == sum_gt(record.gt_from_parent(*v)) ]
			pair = self.select_pair(new_pairs, i, True)
			if pair != (0, 0):
				return pair
			else:
				return self.select_pair(pairs, i, True)
	
	def determine_mat_from(self, i: int) -> int:
		_, mat_gt1, mat_gt2, pat_gt1, pat_gt2 = self.gts(i)
		prev_mat_from = self.from_which_chrom_prev_mat(mat_gt1)
		next_mat_from = self.from_which_chrom_next_mat(mat_gt2)
		# とりあえず、両側Noneはないと仮定する
		if prev_mat_from == 0:
			return next_mat_from
		elif next_mat_from == 0:
			return prev_mat_from
		elif prev_mat_from == next_mat_from:
			return prev_mat_from
		elif self.is_prev_nearer(True):
			return prev_mat_from
		else:
			return next_mat_from
	
	def determine_pat_from(self, i: int) -> int:
		_, mat_gt1, mat_gt2, pat_gt1, pat_gt2 = self.gts(i)
		prev_pat_from = self.from_which_chrom_prev_pat(pat_gt1)
		next_pat_from = self.from_which_chrom_next_pat(pat_gt2)
		# とりあえず、両側Noneはないと仮定する
		if prev_pat_from == 0:
			return next_pat_from
		elif next_pat_from == 0:
			return prev_pat_from
		elif prev_pat_from == next_pat_from:
			return prev_pat_from
		elif self.is_prev_nearer(False):
			return prev_pat_from
		else:
			return next_pat_from
	
	# mat_gt, pat_gt : 0|0 0|1 1|0 1|1を0～3で表す
	def compute_phasing_likelihood(self, mat_gt: int, pat_gt: int) -> float:
		memo = { (0, 0): [0.5, 0.5], (0, 1): [0.9, 0.1], (0, 2): [0.1, 0.9],
				 (1, 0): [0.9, 0.1], (1, 1): [0.99, 0.01], (1, 2): [0.5, 0.5],
				 (2, 0): [0.1, 0.9], (2, 1): [0.5, 0.5], (2, 2): [0.01, 0.99] }
		def probs_from_which_chrom(prev_chrom: int,
								   next_chrom: int) -> list[float]:
			return memo[(prev_chrom, next_chrom)]
		
		def likelihood_each(gt: str, probs_mat: list[float],
									 probs_pat: list[float]) -> float:
			sum_gt = int(gt[0]) + int(gt[2])
			likelihood = sum(probs_mat[i] * probs_pat[j]
								for i, j in product(range(2), repeat=2)
								if ((mat_gt >> i) & 1) +
								   ((pat_gt >> j) & 1) == sum_gt)
			return log(likelihood)
		
		ll = 0.0	# log of likelihood
		for gt, mat_gt1, mat_gt2, pat_gt1, pat_gt2 in self.gen_gts():
			if not Genotype.is_valid(gt, mat_gt, pat_gt):
				ll += log(0.0001)
				continue
			prev_mat_from = self.from_which_chrom_prev_mat(mat_gt1)
			next_mat_from = self.from_which_chrom_next_mat(mat_gt2)
			prev_pat_from = self.from_which_chrom_prev_pat(pat_gt1)
			next_pat_from = self.from_which_chrom_next_pat(pat_gt2)
			probs_mat = probs_from_which_chrom(prev_mat_from, next_mat_from)
			probs_pat = probs_from_which_chrom(prev_pat_from, next_pat_from)
			ll += likelihood_each(gt, probs_mat, probs_pat)
		return ll
	
	def determine_parents_phasing(self):
		if self.record is None:
			return
		
		record = self.record
		lls = []
		for mat_p, pat_p in self.possible_phasings():
			ll = self.compute_phasing_likelihood(mat_p, pat_p)
			lls.append((ll, mat_p, pat_p))
		lls.sort()
		
		mat_gt, pat_gt = self.determine_phasing_core(lls)
		if record is None:
			return
		gt = ['0|0', '1|0', '0|1', '1|1']
		record.v[9] = gt[mat_gt] + record.v[9][3:]
		record.v[10] = gt[pat_gt] + record.v[10][3:]
	
	# どちらから来たか決める
	def select_from(self, froms: list[int], record1: Optional[VCFRecord],
										  record2: Optional[VCFRecord]) -> int:
		if len(froms) == 1:
			return froms[0]
		
		# どちらかは0でない前提
		if froms[0] == 0:
			return froms[1]
		elif froms[1] == 0:
			return froms[0]
		else:	# 両側にRecordがある
			# 近い方を選ぶ
			if self.record is None or record1 is None or record2 is None:
				# ここには来ないはず
				return 0
			elif self.record.pos() * 2 < record1.pos() + record2.pos():
				return froms[0]
			else:
				return froms[1]
	
	def modify_gt(self, i: int) -> str:
		if self.record is None:
			return ""
		
		gt = self.record.get_gt(i)
		prev_mat_from = self.prev_mat_from(i)
		next_mat_from = self.next_mat_from(i)
		prev_pat_from = self.prev_pat_from(i)
		next_pat_from = self.next_pat_from(i)
		if ((prev_mat_from == 0 and next_mat_from == 0) or
								(prev_pat_from == 0 and next_pat_from == 0)):
			return gt
		
		mat_froms = unique_list(prev_mat_from, next_mat_from)
		pat_froms = unique_list(prev_pat_from, next_pat_from)
		pairs = [ x for x in product(mat_froms, pat_froms)
									if x[0] != 0 and x[1] != 0 ]
		mat_from, pat_from = self.select_pair(pairs, i, False)
		if (mat_from, pat_from) == (0, 0):
			if any(mat_from != 0 for mat_from in mat_froms):
				mat_from_nz = self.select_from(mat_froms,
									self.prev_mat_record, self.next_mat_record)
				return self.record.gt_from_mat(mat_from_nz, i+9)
			elif any(pat_from != 0 for pat_from in pat_froms):
				pat_from_nz = self.select_from(pat_froms,
									self.prev_pat_record, self.next_pat_record)
				return self.record.gt_from_pat(pat_from_nz, i+9)
			else:
				return gt	# phasingしない
		else:
			return self.record.gt_from_parent(mat_from, pat_from)
	
	def impute_core(self):
		new_gts = [ self.modify_gt(i)
					for i in range(2, len(self.record.samples)) ]
		self.record.modify_gts(new_gts)
		self.record.modify_parents_type()
	
	def impute(self, necessary_parents_phasing: bool):
		if necessary_parents_phasing:
			self.determine_parents_phasing()
		self.impute_core()
