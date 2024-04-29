from __future__ import annotations

# coding: utf-8
# option.py

from typing import Optional
import sys


#################### Option ####################

class Option:
	def __init__(self, VCF: str, ref: str, ped: str, m: str, r: str,
						families: list[int], chroms: list[int],
						num: int, lp: int, ol: bool, ii: bool,
						ou: bool, ci: bool, out: str):
		self.path_VCF: str						= VCF
		self.path_ped: str						= ped
		self.path_ref_VCF: str					= ref
		self.path_map: str						= m
		self.path_ref: str						= r
		self.families: list[int]				= families
		self.chroms: list[int]					= chroms
		self.num_threads: int					= num
		self.path_out: str						= out
		self.ratio: float						= 0.01
		self.lower_progs: int					= lp
		self.only_large_families: bool			= ol
		self.imputes_isolated_samples: bool		= ii
		self.outputs_unimputed_samples: bool	= ou
		self.corrects_isolated_samples: bool	= ci
	
	def exists_ref(self) -> bool:
		return self.path_ref_VCF != ''
	
	def exists_reference(self) -> bool:
		return self.path_VCF != ''
	
	def print_info(self):
		# required
		print("input VCF : %s" % self.path_VCF, file=sys.stderr)
		print("pedigree : %s" % self.path_ped, file=sys.stderr)
		print("output VCF : %s" % self.path_out, file=sys.stderr)
		
		# optional
		print("number of threads : %s" % self.num_threads, file=sys.stderr)
		print("number of progenies for large family : %s" % self.lower_progs,
															file=sys.stderr)
		if self.path_ref:
			print("ref VCF : %s" % self.path_ref, file=sys.stderr)
		else:
			print("ref VCF is not specified.", file=sys.stderr)
		
		if not self.imputes_isolated_samples:
			print("isolate samples will not be imputed.", file=sys.stderr)
			if self.outputs_unimputed_samples:
				print("but, outputs these samples", file=sys.stderr)
	
	@staticmethod
	def exists(s: str, argv: list[str]) -> bool:
		return s in argv[1:]
	
	@staticmethod
	def flag_value(s: str, argv: list[str]) -> str:
		for option, value in zip(argv, argv[1:]):
			if option == s:
				return value
		else:
			return ''
	
	# str -> [index]
	@staticmethod
	def parse_indices(c: str) -> list[int]:
		indices = []
		v = c.split(',')
		for s in v:
			w = s.split(':')
			if len(w) == 1:
				indices.append(int(w[0]))
			elif len(w) == 2:
				indices.extend(list(range(int(w[0]), int(w[1]))))
			else:
				raise ValueError
		return indices
	
	@staticmethod
	def get_families(argv: list[str]) -> list[int]:
		s = Option.flag_value('-f', argv)
		if s == '':
			return []	# all
		else:
			return Option.parse_indices(s)
	
	@staticmethod
	def get_chroms(argv: list[str]) -> list[int]:
		s = Option.flag_value('-c', argv)
		if s == '':
			return []	# all
		else:
			return Option.parse_indices(s)
	
	@staticmethod
	def get_num_threads(argv: list[str]) -> int:
		s = Option.flag_value('-t', argv)
		if s == '':
			return 1
		else:
			return int(s)
	
	@staticmethod
	def get_lower_progenies(argv: list[str]) -> int:
		s = Option.flag_value('--lower-progs', argv)
		if s == '':
			return 10
		else:
			return int(s)
	
	@staticmethod
	def create(argv: list[str]) -> Optional[Option]:
		try:
			# The top three are required arguments
			path_vcf = Option.flag_value('-i', argv)
			if not path_vcf:
				raise Exception('input VCF not specified.')
			
			path_ped = Option.flag_value('-p', argv)
			if not path_ped:
				raise Exception('pedigree file not specified.')
			
			path_out = Option.flag_value('-o', argv)
			if not path_out:
				raise Exception('output VCF not specified.')
			
			# Optional
			ref_vcf = Option.flag_value('--ref', argv)
			path_map = Option.flag_value('-m', argv)
			path_ref = Option.flag_value('-r', argv)
			families = Option.get_families(argv)
			chroms = Option.get_chroms(argv)
			num_threads = Option.get_num_threads(argv)
			lower_progs = Option.get_lower_progenies(argv)
			ol = Option.exists('-l', argv)
			only_large_families = Option.exists('-large-only', argv)
			impute_isolated = not Option.exists('--not-impute-isolated', argv)
			out_isolated = Option.exists('--out-isolated',  argv)
			if impute_isolated and out_isolated:
				return None
			
			corrects_isolated_samples = Option.exists(
												"--correct-isolated", argv)
			return Option(path_vcf, ref_vcf, path_ped, path_map, path_ref,
							families, chroms, num_threads, lower_progs,
							only_large_families, impute_isolated,
							out_isolated, corrects_isolated_samples, path_out)
		except ValueError:
			return None
		except Exception as e:
			print('error : ' + str(e), file=sys.stderr)
			return None
	
	@staticmethod
	def usage():
		u = ('python graphite.py -i VCF [--ref ref VCF] ' +
				'-p ped [-m map] [-r ref] [-t num_threads] ' +
				'[-f family indices] [-c chrom indices] ' +
				'[--lower-progs lower num progenies] [--large-only] ' +
				'[--not-impute-isolated [--out-isolated]] ' +
				'[--correct-isolated] ' +
				'-o out.')
		messages = ['usage : ' + u,
					'indices: (index|first:last)[,(index|first:last)[,..]]',
					'--large-only: large families only.',
					'--not-impute-isolated: not impute isolated samples.',
					'--out-isolated: output not imputed isolated samples.',
					'--correct-isolated: ' +
						'correct wrong genotypes of isolated samples.']
		for msg in messages:
			print(msg, file=sys.stderr)


#################### OptionImpute ####################

class OptionImpute:
	def __init__(self, max_dist: int, min_positions: int,
									min_graph: int, min_c: float):
		self.max_dist		= max_dist
		self.min_positions	= min_positions
		self.min_graph		= min_graph
		self.min_crossover: float	= 1.0
