from __future__ import annotations

# coding: utf-8

from collections import Counter
import sys

from VCFFamily import *
from VCFImpFamily import FillType, VCFImpFamilyRecord
from TypeDeterminer import ParentComb


#################### VCFHomoHomoRecord ####################

class VCFHomoHomoRecord(VCFImpFamilyRecord):
	def __init__(self, v: list[str], samples: list[str],
						index: int, parents_wrong_type: str, pair: ParentComb):
		super().__init__(v, samples, index, parents_wrong_type, pair)
	
	def is_imputable(self) -> bool:
		return False
	
	def get_fill_type(self) -> FillType:
		if self.pair == ParentComb.P00x11 and self.v[9][:3] == './.':
			return FillType.IMPUTABLE
		else:
			return FillType.FILLED
	
	def gts(self) -> list[str]:
		if self.pair == ParentComb.P00x00:
			return ['0|0'] * len(self.samples)
		elif self.pair == ParentComb.P11x11:
			return ['1|1'] * len(self.samples)
		else:
			if ((self.mat_int_gt() == 0 and self.pat_int_gt() != 0) or
						(self.mat_int_gt() != 2 and self.pat_int_gt() == 2)):
				return ['0|0', '1|1'] + ['0|1'] * self.num_progenies()
			if ((self.mat_int_gt() == 2 and self.pat_int_gt() != 2) or
						(self.mat_int_gt() != 0 and self.pat_int_gt() == 0)):
				return ['1|1', '0|0'] + ['1|0'] * self.num_progenies()
			else:
				return ['./.', './.'] + ['0/1'] * self.num_progenies()
	
	def impute(self) -> VCFHomoHomoRecord:
		v = self.v[:9] + [ gt + s[3:] for gt, s in zip(self.gts(), self.v[9:]) ]
		return VCFHomoHomoRecord(v, self.samples, self.index,
											self.parents_wrong_type, self.pair)


#################### VCFHomoHomo ####################

class VCFHomoHomo(VCFBase, VCFSmallBase, VCFFamilyBase):
	def __init__(self, header: list[list[str]],
						records: list[VCFHomoHomoRecord]):
		VCFBase.__init__(self, header)
		VCFFamilyBase.__init__(self)
		self.records = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def impute(self) -> list[VCFHomoHomo]:
		L: int = self.num_samples()
		records = [ record.impute() for record in self.records ]
		# HeteroHeteroがlistを返すのでそれに合わせる
		return [VCFHomoHomo(self.header, records)]
