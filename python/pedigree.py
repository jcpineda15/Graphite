# coding: utf-8
# pedigree.py

from __future__ import annotations
from typing import Iterator

from common import read_csv


#################### Family ####################

class Family:
	def __init__(self, mat: str, pat: str, progs: list[str]):
		self.mat: str				= mat
		self.pat: str				= pat
		self.progenies: list[str]	= progs
	
	def parents(self) -> tuple[str, str]:
		return (self.mat, self.pat)
	
	def num_progenies(self) -> int:
		return len(self.progenies)
	
	def samples(self) -> list[str]:
		return [self.mat, self.pat] + self.progenies


#################### Progeny ####################

class Progeny:
	def __init__(self, v: list[str]):
		self.family: str	= v[0]
		self.name: str		= v[1]
		self.mat: str		= v[2]
		self.pat: str		= v[3]
	
	def parents(self) -> tuple[str, str]:
		return (self.mat, self.pat)


#################### PedigreeTable ####################

class PedigreeTable:
	def __init__(self, progs: list[Progeny]):
		self.table: list[Progeny]	= progs
	
	def __len__(self) -> int:
		return len(self.table)
	
	def filter_with_parents(self, samples: list[str]) -> PedigreeTable:
		progs = [ p for p in self.table
						if p.mat != '0' and p.pat != '0' and
							p.mat in samples and p.pat in samples ]
		return PedigreeTable(progs)
	
	def generate_parents(self) -> Iterator[tuple[str, str]]:
		return (p.parents() for p in self.table)
	
	def get_children(self, parents: tuple[str, str]) -> list[str]:
		return [ p.name for p in self.table if p.parents() == parents ]
	
	def extract_families(self) -> list[Family]:
		parents = sorted(set((prog.mat, prog.pat) for prog in self.table))
		return [ Family(*p, self.get_children(p)) for p in parents ]
	
	def limit_samples(self, samples: list[str]) -> PedigreeTable:
		set_samples = set(samples)
		progs = [ p for p in self.table if p.name in set_samples ]
		return PedigreeTable(progs)
	
	@staticmethod
	def read(path: str) -> PedigreeTable:
		return PedigreeTable(list(map(Progeny, read_csv(path, ' '))))
