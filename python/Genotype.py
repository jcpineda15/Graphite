from __future__ import annotations

# coding: utf-8
# Genotype.py


#################### Genotype ####################

class Genotype:
	def __init__(self, s: str):
		self.phasing: bool	= s[1] == '|'
		self.gt1: str		= s[0]
		self.gt2: str		= s[2]
	
	def gts(self) -> list[str]:
		return [self.gt1, self.gt2]
	
	def __str__(self) -> str:
		return self.gt1 + ('|' if self.phasing else '/') + self.gt2
	
	@staticmethod
	def int_to_gt(n: int) -> str:
		if n == 0:
			return '0/0'
		elif n == 1:
			return '0/1'
		elif n == 2:
			return '1/1'
		else:
			return './.'
	
	@staticmethod
	def conflicts(mat: Genotype, pat: Genotype, prog: Genotype,
									considers_phasing: bool = True) -> bool:
		if considers_phasing and prog.phasing:
			return not (prog.gt1 in mat.gts() and prog.gt2 in pat.gts())
		else:
			if prog.gt1 == prog.gt2:
				return not (prog.gt1 in mat.gts() and prog.gt2 in pat.gts())
			else:
				return not ((prog.gt1 in mat.gts() and prog.gt2 in pat.gts()) or
							(prog.gt1 in pat.gts() and prog.gt2 in mat.gts()))
	
	@staticmethod
	def is_valid(gt: str, mat_gt: int, pat_gt: int) -> bool:
		mat_gts = Genotype.possible_gts(mat_gt)
		pat_gts = Genotype.possible_gts(pat_gt)
		return len(gt) >= 3 and (gt[0] in mat_gts and gt[2] in pat_gts or
								 gt[2] in mat_gts and gt[0] in pat_gts)
	
	@staticmethod
	def possible_gts(gt: int) -> list[str]:
		if gt == 0:
			return ['0']
		elif gt == 3:
			return ['1']
		else:
			return ['0', '1']
