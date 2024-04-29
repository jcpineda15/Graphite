from __future__ import annotations

# coding: utf-8
# VCFImpFamily.py
# VCFFamilyを継承するこのプロジェクト用のクラス

from abc import ABCMeta, abstractmethod
from functools import reduce
from collections import defaultdict
from operator import add, xor
from enum import Enum

from VCFFamily import *
from TypeDeterminer import ParentComb
from common import is_all_same


#################### FillType ####################

class FillType(Enum):
	MAT = 0
	PAT = 1
	FILLED = 2
	IMPUTABLE = 3
	UNABLE = 4


#################### VCFImpFamilyRecord ####################

class VCFImpFamilyRecord(VCFFamilyRecord, metaclass=ABCMeta):
	def __init__(self, v: list[str], samples: list[str],
						index: int, parents_wrong_type: str, pair: ParentComb):
		super().__init__(v, samples)
		self.index: int = index
		self.parents_wrong_type: str = parents_wrong_type
		self.pair: ParentComb = pair
	
	def is_fixed(self) -> bool:
		return (self.pair in (ParentComb.P00x00, ParentComb.P11x11) or
				(self.pair.is_heterohomo() and self.is_right()))
	
	def is_right(self) -> bool:
		return self.parents_wrong_type == 'Right'
	
	def enable_modification(self):
		self.parents_wrong_type = 'Modifiable'
	
	def set_00x11_parents(self, i: int, gt: int):
		def GT(gt: int) -> str:
			return '{0}|{0}'.format(gt // 2)
		
		GT1 = GT(gt)
		GT2 = GT(2 if gt == 0 else 0)
		j = 1 if i == 0 else 0
		self.set_GT(i, GT1)
		self.set_GT(j, GT2)
		progeny_GT = '0|1' if xor(i == 0, gt == 2) else '1|0'
		for k in range(2, len(self.samples)):
			self.set_GT(k, progeny_GT)
	
	@abstractmethod
	def is_imputable(self) -> bool:
		pass
	
	@abstractmethod
	def get_fill_type(self) -> FillType:
		pass
	
	# あるGenotypeはfixedのみ、その他はfixed以外
	# -> (right GT, [(wrong record, parent index)])
	@staticmethod
	def which_is_fixed(v: list[tuple[VCFImpFamilyRecord, int]]
						) -> tuple[int, list[tuple[VCFImpFamilyRecord, int]]]:
		gts = [ r.get_int_gt(i) for r, i in v ]
		if is_all_same(gts):
			return (gts[0], [])
		
		# collect reads in Genotype
		# ここで./.を無視してもよいはず
		dic = defaultdict(list)	# { int gt: [(record, parent index)] }
		for (r, i), gt in zip(v, gts):
			dic[gt].append((r, i))
		items = list(dic.items())	# [(int gt, [(record, parent index)])]
		
		bs = [ all(r.is_fixed() for r, i in w) for gt, w in items ]
		# 全部fixedが一つでないとあきらめる
		# 先頭のGenotypeを返すが、空のレコードのリストを返すので何でもよい
		if sum(1 for b in bs if b) != 1:
			return (gts[0], [])
		
		bs2 = [ all(not r.is_fixed() for r, i in w) for gt, w in items ]
		if sum(1 for b in bs2 if not b) != 1:
			return (gts[0], [])
		
		fixed_index = next(i for i, b in enumerate(bs) if b)
		fixed_GT = items[fixed_index][0]
		if fixed_GT == 1:	# 0/1
			return (gts[0], [])
		
		# 0/0 x 1/1パターン以外のRecordはあとで直す機会がある
		wrongs = [ (r, k) for j, item in enumerate(items) if j != fixed_index
							for r, k in item[1] if r.pair == ParentComb.P00x11 ]
		return (fixed_GT, wrongs)
	
	# 子どものGenotypeから0/0 x 1/1と推定される場合、
	# 親のどちらが0/0でどちらが1/1なのかを間違えている場合が多い
	# 親が複数の家系にまたがっているとき、修正できるかもしれない
	@staticmethod
	def modify_00x11(records: list[VCFImpFamilyRecord]):
		# { sample: [(VCFImpFamilyRecord, mat or pat)] }
		dic = defaultdict(list)
		for r in records:
			dic[r.samples[0]].append((r, 0))
			dic[r.samples[1]].append((r, 1))
		
		for v in dic.values():
			if len(v) < 2:
				continue
			
			right_gt, wrongs = VCFImpFamilyRecord.which_is_fixed(v)
			for record, i in wrongs:
				record.set_00x11_parents(i, right_gt)
