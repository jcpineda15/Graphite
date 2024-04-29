from __future__ import annotations

# coding: utf-8
# SampleManager.py

from pedigree import Family


#################### KnownFamily ####################

# 親がsampleにあるかどうかを記憶している
class KnownFamily(Family):
	def __init__(self, mat: str, pat: str,
					mat_known: bool, pat_known: bool, progs: list[str]):
		super().__init__(mat, pat, progs)
		self.mat_known: bool = mat_known
		self.pat_known: bool = pat_known
	
	def is_both_unknown(self) -> bool:
		return not self.mat_known and not self.pat_known
	
	def is_one_unknown(self) -> bool:
		return ((not self.mat_known and self.pat_known) or
				(self.mat_known and not self.pat_known))
	
	def known_parents(self) -> list[str]:
		parents = []
		if self.mat_known:
			parents.append(self.mat)
		if self.pat_known:
			parents.append(self.pat)
		return parents
