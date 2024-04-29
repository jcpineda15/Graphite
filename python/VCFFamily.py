from __future__ import annotations

# coding: utf-8
# VCFFamily.py

from abc import ABC, abstractmethod
from itertools import count, groupby

from typing import Dict, Iterator, Optional, IO, TypeVar
from VCF import *
from pedigree import *


#################### VCFFamilyRecord ####################

class VCFFamilyRecord(VCFRecord):
	def __init__(self, v: list[str], samples: list[str]):
		super().__init__(v, samples)
	
	def mat_gt(self) -> str:
		return self.v[9]
	
	def pat_gt(self) -> str:
		return self.v[10]
	
	def mat_GT(self) -> str:
		return self.get_GT(0)
	
	def pat_GT(self) -> str:
		return self.get_GT(1)
	
	def mat_int_gt(self) -> int:
		return self.get_int_gt(0)
	
	def pat_int_gt(self) -> int:
		return self.get_int_gt(1)
	
	def progeny_gts(self) -> list[int]:
		return [ self.get_int_gt(i) for i in range(2, len(self.v) - 9) ]
	
	def progeny_int_gts(self) -> list[int]:
		return [ self.get_int_gt(i) for i in range(2, len(self.v) - 9) ]
	
	def num_progenies(self) -> int:
		return len(self.v) - 11


#################### VCFFamilyBase ####################

# Inherite this not VCFFamily for A family VCF
class VCFFamilyBase(ABC):
	def __init__(self):
		pass
	
	@abstractmethod
	def get_header(self) -> list[list[str]]:
		pass
	
	@abstractmethod
	def get_samples(self) -> list[str]:
		pass
	
	@abstractmethod
	def __len__(self) -> int:
		pass
	
	@abstractmethod
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		pass
	
	def mat(self) -> str:
		return self.get_samples()[0]
	
	def pat(self) -> str:
		return self.get_samples()[1]
	
	def parents(self) -> tuple[str, str]:
		return (self.mat(), self.pat())
	
	def num_progenies(self) -> int:
		return len(self.get_samples()) - 2


#################### VCFFamily ####################

class VCFFamily(VCFBase, VCFSmallBase, VCFFamilyBase):
	def __init__(self, header: list[list[str]], records: list[VCFFamilyRecord]):
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		VCFFamilyBase.__init__(self)
		self.records = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	@staticmethod
	def create(vcf: VCFSmall, samples: list[str]):
		columns = vcf.extract_columns(samples)
		h = vcf.header[-1]
		header = vcf.header[:-1] + [h[:9] + [ h[c] for c in columns ]]
		# [VCFFamilyRecord]
		records = [ VCFFamily.subset(r, columns) for r in vcf.records ]
		return VCFFamily(header, records)
	
	# vcf1にあるsampleはvcf1から、
	# そうでないsampleはvcf2からGenotypeを取って新たなVCFを作る
	@staticmethod
	def create_by_two_vcfs(vcf1: VCFSmallBase, vcf2: VCFSmallBase,
										samples: list[str]) -> VCFFamily:
		# samplesは[mat, pat, prog1, ...]の前提
		columns1 = vcf1.extract_columns(samples)
		columns2 = vcf2.extract_columns(samples)
		new_header = vcf1.trim_header(samples)
		new_records = []
		for i in range(len(vcf1)):
			record1 = vcf1.get_record(i)
			record2 = vcf2.get_record(i)
			v = record1.v[:9]
			for i in range(len(samples)):
				if columns1[i] != -1:
					v.append(record1.v[columns1[i]])
				elif columns2[i] != -1:
					v.append(record2.v[columns2[i]])
				else:	# どちらのVCFにもそのサンプルは無い
					v.append('./.')
			new_record = VCFFamilyRecord(v, samples)
			new_records.append(new_record)
		return VCFFamily(new_header, new_records)
	
	@staticmethod
	def subset(record: VCFRecord, columns: list[int]) -> VCFFamilyRecord:
		v = record.v
		new_v = v[:9] + [ v[c] if c != -1 else './.' for c in columns ]
		samples = [ record.samples[c-9] if c != -1 else '0'
												for c in columns ]
		return VCFFamilyRecord(new_v, samples)
	
	# [VCFFamily] -> VCFFamily
	@staticmethod
	def join(vcfs):
		records = [ record for vcf in vcfs for record in vcf.records ]
		return VCFFamily(vcfs[0].header, records)
	
	@staticmethod
	def convert(vcf: VCFFamilyBase) -> VCFFamily:
		records: list[VCFFamilyRecord] = [ vcf.get_family_record(i)
													for i in range(len(vcf)) ]
		return VCFFamily(vcf.get_header(), records)

