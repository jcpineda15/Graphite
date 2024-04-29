
# coding: utf-8
# VCFFillableRecord.py

from __future__ import annotations
from functools import reduce
from itertools import *
from collections import Counter

from VCFFamily import *
from VCFImpFamily import FillType, VCFImpFamilyRecord
from TypeDeterminer import ParentComb
from common import *


#################### VCFFillableRecord ####################

class VCFFillableRecord(VCFFamilyRecord):
	def __init__(self, v: list[str], samples: list[str], index: int,
											type: FillType, pair: ParentComb):
		super().__init__(v, samples)
		self.index: int = index
		self.type: FillType = type
		self.pair: ParentComb = pair
	
	def group_ids(self) -> list[int]:
		w = self.v[7].split('=')
		return list(map(int, w[1].split(',')))
	
	def set_group_ids(self, ids: list[int]):
		self.v[7] = 'GR=%s' % (','.join(map(str, ids)))
	
	def gt_from_parent(self, mat_from: int, pat_from: int) -> str:
		gt_from_mat = '.' if mat_from == 0 else self.v[9][mat_from*2-2]
		gt_from_pat = '.' if pat_from == 0 else self.v[10][pat_from*2-2]
		return gt_from_mat + '|' + gt_from_pat
	
	def gt_from_mat(self, mat_from: int, c: int) -> str:
		gt = self.get_GT(c - 9)
		int_gt = self.get_int_gt(c - 9)
		mat_gt = self.v[9][mat_from*2-2]
		pat_GT = self.pat_GT()
		if int_gt == 0:
			if mat_gt == '0':
				return '0|0'
			elif '0' in pat_GT:
				return '1|0'
			else:
				return '1|1'
		elif int_gt == 2:
			if mat_gt == '1':
				return '1|1'
			elif '1' in pat_GT:
				return '0|1'
			else:
				return '0|0'
		elif int_gt == 1:
			if mat_gt == '0':
				return '0|1'
			else:
				return '1|0'
		else:
			return '.|.'
	
	def gt_from_pat(self, pat_from: int, c: int) -> str:
		gt = self.get_GT(c - 9)
		int_gt = self.get_int_gt(c - 9)
		pat_gt = self.v[10][pat_from*2-2]
		mat_GT = self.mat_GT()
		if int_gt == 0:
			if pat_gt == '0':
				return '0|0'
			elif '0' in mat_GT:
				return '0|1'
			else:
				return '1|1'
		elif int_gt == 2:
			if pat_gt == '1':
				return '1|1'
			elif '1' in mat_GT:
				return '1|0'
			else:
				return '0|0'
		elif int_gt == 1:
			if pat_gt == '0':
				return '1|0'
			else:
				return '0|1'
		else:
			return '.|.'
	
	def is_mat_hetero(self) -> bool:
		return self.v[9][0] != self.v[9][2]
	
	def mat_from(self, c: int) -> int:
		if self.v[c][0] == '.':
			return 0
		elif not self.is_hetero(0):
			return 0
		elif self.v[c][0] == self.v[9][0]:
			return 1
		else:
			return 2
	
	def pat_from(self, c: int) -> int:
		if self.v[c][0] == '.':
			return 0
		elif not self.is_hetero(1):
			return 0
		elif self.v[c][2] == self.v[10][0]:
			return 1
		else:
			return 2
	
	def find_geno_type(self, type: str) -> int:
		for i, t in enumerate(self.v[8].split(':')):
			if t == type:
				return i
		else:
			return -1
	
	def fill_PGT(self) -> None:
		i_GT: int = self.find_geno_type('GT')
		assert(i_GT != -1)
		i_PGT: int = self.find_geno_type('PGT')
		if i_PGT == -1:
			return
		
		for j in range(9, len(self.v)):
			v = self.v[j].split(':')
			if i_PGT >= len(v):		# 補完の際にGTだけになった
				continue
			v[i_PGT] = v[i_GT]
			self.v[j] = ':'.join(v)
	
	def inverse_gt(self, gt: str, inv: bool) -> str:
		return ("1" if gt == "0" else "0") if inv else gt
	
	def inverse_prog_gt(self, gt: str, inv_mat: bool, inv_pat: bool) -> str:
		return (self.inverse_gt(gt[0], inv_mat) + "|" +
				self.inverse_gt(gt[2], inv_pat) + gt[3:])
	
	def inverse_prog_gts(self, prog_gts: list[str],
							inv_mat: bool, inv_pat: bool) -> list[str]:
		inv_prog_gts = [ self.inverse_prog_gt(gt, inv_mat, inv_pat)
												for gt in prog_gts ]
		return inv_prog_gts
	
	def inverse_parents_gts(self, inv_mat: bool, inv_pat: bool):
		# both must be hetero
		if inv_mat:
			self.set_GT(0, self.v[9][2] + "|" + self.v[9][0])
		if inv_pat:
			self.set_GT(1, self.v[10][2] + "|" + self.v[10][0])
	
	def is_same_gts(self, gt1: str, gt2: str) -> bool:
		if gt2 == "0/1":
			return gt1 == "0|1" or gt1 == "1|0"
		elif gt2 == "0/0":
			return gt1 == "0|0"
		elif gt2 == "1/1":
			return gt1 == "1|1"
		else:
			return False;
	
	def is_near_prog_gts(self, gts: list[str]) -> bool:
		num = 0
		dist = 0
		for i in range(len(gts)):
			if self.v[i+11] != "0/1":
				num += 1
			if not self.is_same_gts(gts[i], self.v[i+11]):
				dist += 1
		return dist < max(1, num // 2)
	
	def modify_gts(self, new_prog_gts: list[str]):
		for inv_mat, inv_pat in product((False, True), repeat=2):
			inv_prog_gts = self.inverse_prog_gts(new_prog_gts, inv_mat, inv_pat)
			if self.is_near_prog_gts(inv_prog_gts):
				self.inverse_parents_gts(inv_mat, inv_pat)
				for c in range(11, len(self.v)):
					self.v[c] = inv_prog_gts[c-11]
				return
		for c in range(11, len(self.v)):
			self.v[c] = new_prog_gts[c-11]

	def modify_parents_type(self):
		if (self.pair != ParentComb.P00x11 and
				((self.v[9][:3] == '0|0' and self.v[10][:3] == '1|1') or
				 (self.v[9][:3] == '1|1' and self.v[10][:3] == '0|0'))):
			self.pair = ParentComb.P00x11
	
	def from_which_chrom(self, i: int, mat: bool) -> int:
		j = 0 if mat else 1
		parent_gt = self.v[j+9]
		return 1 if parent_gt[0] == self.get_gt(i)[j*2] else 2
	
	@staticmethod
	def convert(record: VCFImpFamilyRecord) -> VCFFillableRecord:
		type = record.get_fill_type()
		return VCFFillableRecord(record.v, record.samples, record.index,
															type, record.pair)
	
	@staticmethod
	def merge(records: list[VCFFillableRecord],
								samples: list[str]) -> VCFRecord:
		v = records[0].v[:]
		for record in records[1:]:
			v.extend(record.v[9:])
		return VCFRecord(v, samples)
	
	def decide_by_majority(self, GTs: list[str]) -> str:
		# 多数決できるなら多数決
		c = Counter(GTs)
		dic_num = Counter(c.values())
		max_num = max(dic_num.keys())
		if dic_num[max_num] == 1:	# 最多のGenotypeは一つ
			for GT, num in c.items():
				if num == max_num:
					return GT
			else:
				# ここには来ないはず
				return GT
		else:						# 複数ある場合は乱数的なものを使う
			candidate_GTs = []
			for GT, num in c.items():
				if num == max_num:
					candidate_GTs.append(GT)
			candidate_GTs.sort()
			# posをhash化してGTを選ぶ
			d = len(candidate_GTs)
			i = reduce(lambda x, y: x + int(y), self.v[1], 0) % d
			return candidate_GTs[i]
	
	def swap_parents(self, pos: int, GT: str):
		if GT == self.v[pos+9][:3]:
			return
		elif GT in ('0|0', '1|1'):
			is_mat_00 = (pos == 0) ^ (GT == '1|1')
			self.set_GT(0, '0|0' if is_mat_00 else '1|1')
			self.set_GT(1, '1|1' if is_mat_00 else '0|0')
			# 子どもも入れ換えなければならない
			prog_GT = '0|1' if is_mat_00 else '1|0'
			for i in range(2, len(self.samples)):
				self.set_GT(i, prog_GT)
	
	@staticmethod
	def decide_duplicated_Genotype(records: list[VCFFillableRecord],
											ps: list[tuple[int, int]]) -> str:
		# 全部./.ならそのまま
		# 子どもがいたらそれ
		# ./.を除いて全部同じだったらそれ
		# ./.と0/0 x 1/1で全部なら0/0 x 1/1で多数決
		# ./.と0/0 x 1/1を除いて全部同じだったらそれ
		# ./.と0/0 x 1/1を除いて複数種あれば多数決
		GTs = [ records[i].v[j+9][:3] for i, j in ps ]
		if all(GT == './.' for GT in GTs):
			return './.'
		
		# 子どもが優先
		for (i, j), GT in zip(ps, GTs):
			if j >= 2 and GT != './.':
				return GT
		
		GTs2 = [ GT for GT in GTs if GT != './.' ]
		if is_all_same(GTs2):
			return GTs2[0]
		
		GTs3 = [ GT for (i, j), GT in zip(ps, GTs)
					if GT != './.' and records[i].pair != ParentComb.P00x11 ]
		if not GTs3:	# 0/0 x 1/1しかない
			return records[0].decide_by_majority(GTs2)
		elif is_all_same(GTs3):
			return GTs3[0]
		else:
			return records[0].decide_by_majority(GTs3)
	
	@staticmethod
	def integrate_each_sample(records: list[VCFFillableRecord],
												ps: list[tuple[int, int]]):
		GT = VCFFillableRecord.decide_duplicated_Genotype(records, ps)
		
		for i, j in ps:
			# 0/0 x 1/1の親なら相手と入れ替えなければならない
			if records[i].pair == ParentComb.P00x11 and j <= 1:
				records[i].swap_parents(j, GT)
	
	@staticmethod
	def integrate(records: list[VCFFillableRecord], samples: list[str],
						pos_samples: list[list[tuple[int, int]]]) -> VCFRecord:
		# 0|0 x 1|1ならひっくり返せる
		for ps in pos_samples:
			GTs = [ records[i].v[j+9][:3] for i, j in ps ]
			if not is_all_same(GTs):
				VCFFillableRecord.integrate_each_sample(records, ps)
		
		# ひっくり返した後にGTをまとめる
		v: list[str] = records[0].v[:9]
		for ps in pos_samples:
			i, j = ps[0]
			v.append(records[i].v[j+9])
		return VCFRecord(v, samples)
