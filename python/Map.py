from __future__ import annotations

# coding: utf-8
# Map.py

from itertools import *
from typing import Iterator, List, Any, IO

from common import *


#################### Map ####################

class Map:
	class MapRecord:
		def __init__(self, v: list[Any]):
			self.chr: str	= v[0]
			self.cM: float	= float(v[1])
			self.Mbp: float	= float(v[2])
	
	def __init__(self, records: List['MapRecord']):
		self.records: List['Map.MapRecord'] = records
		self.chr_maps: List[Map] = []
	
	def is_empty(self) -> bool:
		return not self.records
	
	def total_cM(self) -> float:
		return self.records[-1].cM
	
	def bp_to_cM(self, bp: float) -> float:
		Mbp: float = bp * 1e-6
		if self.is_empty():
			return Mbp
		
		index: int = self.binary_search(0, len(self.records)-1, Mbp)
		r1: Map.MapRecord = self.records[index]
		r2: Map.MapRecord = self.records[index+1]
		if r1.Mbp == r2.Mbp:
			return r1.cM
		return r1.cM + (r2.cM-r1.cM) / (r2.Mbp-r1.Mbp) * (Mbp-r1.Mbp)
	
	def binary_search(self, first: int, last: int, Mbp: float) -> int:
		if first == last - 1:
			return first
		
		mid: int = (first + last) // 2
		if Mbp < self.records[mid].Mbp:
			return self.binary_search(first, mid, Mbp)
		else:
			return self.binary_search(mid, last, Mbp)
	
	def create_chr_maps(self) -> list[Map]:
		if self.is_empty():
			# Make one map and use it for every chromosome.
			# That makes it easier to delete.
			return []
		else:
			return list(self.divide_into_chromosomes())
	
	def iter_chr_maps(self) -> Iterator[Map]:
		if self.is_empty():
			while True:
				yield self
		else:
			for m in self.chr_maps:
				yield m
	
	def divide_into_chromosomes(self) -> Iterator[Map]:
		for chr, v in groupby(self.records, key=lambda r: r.chr):
			yield Map(list(v))
	
	def display_info(self, out: IO):
		print("Genetic Map : ", end='', file=out)
		if self.is_empty():
			print("default map(1Mbp=1cM).", file=out)
		else:
			args = (len(self.chr_maps), self.total_cM())
			print("%d chrmosomes %d cM." % args, file=out)
	
	@staticmethod
	def read(path: str) -> Map:
		if path:
			records = [ Map.MapRecord(v) for v in read_csv(path) ]
		else:
			records = []
		m = Map(records)
		m.chr_maps = m.create_chr_maps()
		return m
	
	@staticmethod
	def default_map() -> Map:
		r1 = Map.MapRecord(['1', 0.0, 0.0])
		r2 = Map.MapRecord(['1', 1.0, 1.0])
		return Map([r1, r2])


#################### VCFMeasurable ####################

class VCFMeasurable(object):
	def __init__(self, map_: Map):
		self.map: Map = map_
	
	def cM(self, pos: int) -> float:
		return self.map.bp_to_cM(pos)
