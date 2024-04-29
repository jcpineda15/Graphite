from __future__ import annotations

# coding: utf-8
# imputeしないバージョン
# selfのときのみimputeするバージョンを使う

from VCFFamily import *
from VCFImpFamily import FillType, VCFImpFamilyRecord
from TypeDeterminer import ParentComb


#################### VCFHeteroHeteroLiteRecord ####################

class VCFHeteroHeteroLiteRecord(VCFImpFamilyRecord):
	def __init__(self, v: list[str], samples: list[str],
								index: int, parents_wrong_type: str):
		super().__init__(v, samples, index,
								parents_wrong_type, ParentComb.P01x01)
	
	def is_imputable(self) -> bool:
		return False
	
	def get_fill_type(self) -> FillType:
		return FillType.IMPUTABLE


#################### VCFHeteroHeteroLite ####################

class VCFHeteroHeteroLite(VCFBase, VCFFamilyBase):
	def __init__(self, header: list[list[str]],
						records: list[VCFHeteroHeteroLiteRecord]):
		VCFBase.__init__(self, header)
		VCFFamilyBase.__init__(self)
		self.records = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
