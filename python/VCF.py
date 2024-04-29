# coding: utf-8
# VCF.py

from __future__ import annotations
from typing import Dict, Generator, Iterator, Optional, IO, TypeVar
from abc import ABC, abstractmethod
from itertools import *
from math import log
import re
import gzip

from common import *


#################### library ####################


#################### VCFRecord ####################

class VCFRecord(object):
	def __init__(self, v: list[str], samples: list[str]):
		self.v: list[str] = v
		self.samples: list[str] = samples
	
	def copy(self) -> VCFRecord:
		return VCFRecord(self.v[:], self.samples)
	
	def position(self) -> tuple[str, int]:
		return (self.v[0], int(self.v[1]))
	
	def chrom(self) -> str:
		return self.v[0]
	
	def pos(self) -> int:
		return int(self.v[1])
	
	def format(self) -> str:
		return self.v[8]
	
	def gt(self, sample: list[str]) -> Optional[str]:
		for gt, sample_ in zip(self.v[9:], self.samples):
			if sample_ == sample:
				return gt
		else:
			return None
	
	def is_NA(self, i):
		return self.v[i+9][0] == '.' or self.v[i+9][2] == '.'
	
	def gts(self) -> list[str]:
		return self.v[9:]
	
	def get_GT(self, i: int) -> str:
		return self.v[i+9][:3]
	
	def get_int_gt(self, i: int) -> int:
		s = self.v[i+9]
		try:
			return int(s[0]) + int(s[2])
		except ValueError:
			return -1
	
	def get_int_gts(self) -> list[int]:
		return [ self.get_int_gt(i) for i in range(len(self.samples)) ]
	
	def is_homo(self, i: int) -> bool:
		s = self.v[9+i]
		return s[0] == s[2]
	
	def is_hetero(self, i: int) -> bool:
		gt = self.v[i+9]
		return gt[0] != '.' and gt[0] != gt[2]
	
	def get_gt(self, i: int) -> str:
		return self.v[i+9]
	
	def write(self, out: IO):
		write_tsv(self.v, out)
	
	def set_GT(self, i: int, GT: str):
		self.v[i+9] = GT + self.v[i+9][3:]


#################### VCFBase ####################

class VCFBase(object):
	def __init__(self, header: list[list[str]]):
		self.header: list[list[str]] = header
		self.samples: list[str] = self.header[-1][9:]
		self.dic: Dict[str, int] = dict(zip(self.samples, count()))
		self.chrs: Dict[str, int] = { }
		self.__determine_chromosome_id()
	
	def __determine_chromosome_id(self):
		pat = re.compile(r'##contig=<ID=(.+),length=(\d+)>')
		id = 0
		for s in (v[0] for v in self.header):
			m = pat.match(s)
			if m:
				self.chrs[m.group(1)] = id
				id += 1
	
	def write_header(self, out: IO):
		for v in self.header:
			write_tsv(v, out)
	
	def get_header(self) -> list[list[str]]:
		return self.header
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def record_position(self, record: VCFRecord) -> tuple[int, int]:
		return self.position(record.position())
	
	def position(self, position) -> tuple[int, int]:
		chr, pos = position
		if chr not in self.chrs:
			self.chrs[chr] = len(self.chrs) + 1
		return (self.chrs[chr], pos)
	
	def chr(self, chr_id: int) -> str:
		for chr, id in self.chrs.items():
			if id == chr_id:
				return chr
		return ''
	
	def num_samples(self) -> int:
		return len(self.samples)
	
	@staticmethod
	def read_header(g: Iterator[list[str]]) -> Generator[list[str], None, None]:
		for v in g:
			yield v
			if v[0].startswith('#CHROM'):
				break


#################### VCFHuge ####################

# 巨大ファイルなので、1行ずつ読んでは捨てる前提
class VCFHuge(VCFBase):
	def __init__(self, header: list[list[str]], g: Iterator[list[str]]):
		super().__init__(header)
		self.g: Iterator[list[str]] = g
	
	def __iter__(self):
		return self
	
	def __next__(self):
		record = VCFRecord(next(self.g), self.samples)
		self.record_position(record)	# self.chrsを作るために必要
		return record
	
	# posまで進める
	def proceed(self, pos: tuple[int, int]) -> VCFRecord:
		while True:
			record = next(self)
			pos2 = self.record_position(record)
			if pos2 == pos:
				return record
	
	def close(self):
		self.g.close()
	
	def get_header(self):
		return self.header
	
	def divide_into_chromosomes(self) -> Generator[VCFSmall, None, None]:
		for chr, v in groupby(self, key=VCFRecord.chrom):
			yield VCFSmall(self.get_header(), list(v))
	
	@staticmethod
	def read(path) -> VCFHuge:
		gen_vecs = read_tsv(path)
		header = list(VCFBase.read_header(gen_vecs))
		return VCFHuge(header, gen_vecs)
	
	@staticmethod
	def simple_header(samples: list[str]) -> list[str]:
		return ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
						"FILTER", "INFO", "FORMAT"] + samples


#################### VCFSmallBase ####################

class VCFSmallBase(ABC):
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
	def get_record(self, i: int) -> VCFRecord:
		pass
	
	def trim_header(self, samples: list[str]) -> list[list[str]]:
		header = self.get_header()
		return header[:-1] + [header[-1][:9] + samples]
	
	def extract_columns(self, samples: list[str]) -> list[int]:
		dic = dict(zip(self.get_samples(), count(9)))
		return [ dic.get(sample, -1) for sample in samples ]
	
	def extract_samples(self, samples) -> VCFSmall:
		header = self.trim_header(samples)
		cs = self.extract_columns(samples)
		new_records: list[VCFRecord] = []
		for i in range(len(self)):
			record = self.get_record(i)
			v = record.v
			new_v = v[:9] + [ v[c] for c in cs ]
			new_record = VCFRecord(new_v, samples)
			new_records.append(new_record)
		return VCFSmall(header, new_records)
	
	def clip_raw_haplotype(self, sample_index: int, side: int) -> list[int]:
		hap: list[int] = []
		for j in range(len(self)):
			record = self.get_record(j)
			gt = record.get_GT(sample_index)
			c = gt[side*2]
			hap.append(-1 if c == '.' else int(c))
		return hap
	
	def write_header(self, out: IO):
		for v in self.get_header():
			write_tsv(v, out)
	
	def write(self, out: IO, with_header: bool=True):
		if with_header:
			self.write_header(out)
		for i in range(len(self)):
			self.get_record(i).write(out)


#################### VCFSmall ####################

class VCFSmall(VCFBase, VCFSmallBase):
	def __init__(self, header: list[list[str]], records: list[VCFRecord]):
		super().__init__(header)
		self.records: list[VCFRecord] = records
		for record in self.records:
			if record is not None:
				self.record_position(record)
	
	def __len__(self) -> int:
		return len(self.records)
	
	def copy(self) -> VCFSmall:
		records = [ r.copy() for r in self.records ]
		return VCFSmall(self.header, records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def divide_into_chromosomes(self) -> Iterator[VCFSmall]:
		for chr, v in groupby(self.records, key=VCFRecord.chrom):
			yield VCFSmall(self.header, list(v))
	
	@staticmethod
	def read(path: str) -> VCFSmall:
		g_vecs = read_tsv(path)
		header = list(VCFBase.read_header(g_vecs))
		samples = header[-1][9:]
		records = [ VCFRecord(v, samples) for v in g_vecs ]
		return VCFSmall(header, records)
	
	# samplesの順番で結合
	@staticmethod
	def join(vcfs: list[VCFSmallBase], samples: list[str]) -> VCFSmall:
		dic = { }
		for vcf in vcfs:
			# sに'0'があっても問題ない
			for c, s in enumerate(vcf.get_samples(), 9):
				dic[s] = (vcf, c)
		cols = []
		for s in samples:
			if s in dic:
				cols.append((s,) + dic[s])
		
		new_samples = [ s for s, vcf, c in cols ]
		new_records = []
		for i in range(len(vcfs[0])):
			v = vcfs[0].get_record(i).v[:9]
			for s, vcf, c in cols:
				v.append(vcf.get_record(i).v[c])
			new_records.append(VCFRecord(v, new_samples))
		new_header = vcfs[0].trim_header(new_samples)
		return VCFSmall(new_header, new_records)
	
	@staticmethod
	def convert(vcf: VCFSmallBase) -> VCFSmall:
		records: list[VCFRecord] = [ vcf.get_record(i)
											for i in range(len(vcf)) ]
		return VCFSmall(vcf.get_header(), records)


#################### main ####################

__all__ = ['VCFRecord', 'VCFBase', 'VCFHuge', 'VCFSmallBase', 'VCFSmall']
