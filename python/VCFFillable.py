from __future__ import annotations

# coding: utf-8
# VCFFillable.py

from itertools import *
from collections import defaultdict

from typing import Optional

from VCFFamily import *
from VCFImpFamily import FillType, VCFImpFamilyRecord
from VCFFillableRecord import VCFFillableRecord
from VCFHeteroHomo import *
from group import Groups
from RecordSet import RecordSet
from Genotype import Genotype
from option import *
from common import *


#################### VCFFillable ####################

class VCFFillable(VCFBase, VCFSmallBase, VCFFamilyBase):
	def __init__(self, header: list[list[str]],
							records: list[VCFFillableRecord]):
		VCFBase.__init__(self, header)
		VCFFamilyBase.__init__(self)
		self.records: list[VCFFillableRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def modify(self) -> None:
		# FillTypeでrecordを分ける
		groups = Groups.create(self.records)
		for record_set in groups.generate_record_sets():
			record_set.determine_parents_phasing()
			self.__impute_core(record_set)
		
		for i, record in enumerate(self.records):
			# 家系ごとで./.にしたGenotypeを補完
			if record.type == FillType.MAT:
				self.__impute_NA_mat(i)
			elif record.type == FillType.PAT:
				self.__impute_NA_pat(i)
			elif record.type in (FillType.IMPUTABLE, FillType.UNABLE):
				self.__impute_others(i)
		
		for record in self.records:
			record.fill_PGT()
	
	def phase_hetero_hetero(self) -> None:
		# FillTypeでrecordを分ける
		groups = Groups.create(self.records)
		for record_set in groups.generate_record_sets():
			record_set.impute(False)
#			self.__impute_core(record_set)
		
		# この部分はフルで必要？
		for i, record in enumerate(self.records):
			# 家系ごとで./.にしたGenotypeを補完
			if record.type == FillType.MAT:
				self.__impute_NA_mat(i)
			elif record.type == FillType.PAT:
				self.__impute_NA_pat(i)
			elif record.type in (FillType.IMPUTABLE, FillType.UNABLE):
				self.__impute_others(i)
	
	def __impute_core(self, record_set: RecordSet):
		def gts(record: Optional[VCFFillableRecord]) -> list[str]:
			return record.v[11:] if record else [''] * (len(self.samples) - 2)
		
		def sum_gt(gt: str) -> int:
			try:
				return int(gt[0]) + int(gt[2])
			except ValueError:
				return -1
		
		def both_nearests() -> tuple[int, int]:
			# ここに来るときは、全てNoneではない
			if (record is None or mat_record1 is None or mat_record2 is None
							   or pat_record1 is None or pat_record2 is None):
				return (0, 0)
			if record_set.is_prev_nearer(True):
				mat_from = prev_mat_from
			else:
				mat_from = next_mat_from
			if record.pos() * 2 < pat_record1.pos() + pat_record2.pos():
				pat_from = prev_pat_from
			else:
				pat_from = next_pat_from
			return (mat_from, pat_from)
		
		# [(mat_from, pat_from)] -> (mat_from, pat_from)
		def nearest_froms(pairs: list[tuple[int, int]]) -> tuple[int, int]:
			if len(pairs) == 4:		# N/Aのときにまれにあり得る
				return both_nearests()
			if pairs[0][0] == pairs[1][0]:		# matが同じ
				if record is None or pat_record1 is None or pat_record2 is None:
					return (0, 0)
				if record.pos() * 2 < pat_record1.pos() + pat_record2.pos():
					return (pairs[0][0], prev_pat_from)
				else:
					return (pairs[0][0], next_pat_from)
			elif pairs[0][1] == pairs[1][1]:	# patが同じ
				if record is None or mat_record1 is None or mat_record2 is None:
					return (0, 0)
				if record.pos() * 2 < mat_record1.pos() + mat_record2.pos():
					return (prev_mat_from, pairs[0][1])
				else:
					return (next_mat_from, pairs[0][1])
			else:	# 両親とも乗り換えている（滅多にない）
				return both_nearests()
		
		def select_pair(pairs: list[tuple[int, int]], gt: str,
									selected: bool = False) -> tuple[int, int]:
			if record is None:
				return (0, 0)
			elif not pairs:
				return (0, 0)
			elif len(pairs) == 1:
				return pairs[0]
			elif gt[0] == '.':
				return nearest_froms(pairs)
			elif selected:
				return nearest_froms(pairs)
			else:
				new_pairs = [ v for v in pairs
							if sum_gt(gt) == sum_gt(record.gt_from_parent(*v)) ]
				pair = select_pair(new_pairs, gt, True)
				if pair != (0, 0):
					return pair
				else:
					return select_pair(pairs, gt, True)
		
		def is_same_gts(gt1: str, gt2: str) -> bool:
			if gt2 == '0/1':
				return gt1 == '0|1' or gt1 == '1|0'
			elif gt2 == '0/0':
				return gt1 == '0|0'
			elif gt2 == '1/1':
				return gt1 == '1|1'
			else:
				return False
		
		def is_near(gts1: list[str], gts2: list[str]) -> bool:
			num = sum(1 for gt2 in gts2 if gt2 != '0/1')
			dist = sum(1 for gt1, gt2 in zip(gts1, gts2)
										if not is_same_gts(gt1, gt2))
			return dist < max(1, num // 2)
		
		def inverse_gts(gts: list[str], inv_mat: bool,
										inv_pat: bool) -> list[str]:
			def inverse(s: str) -> str:
				return '1' if s == '0' else '0'
			
			if inv_mat:
				if inv_pat:
					return [ inverse(gt[0]) + '|' + inverse(gt[2]) + gt[3:]
																for gt in gts ]
				else:
					return [ inverse(gt[0]) + '|' +  gt[2:] for gt in gts ]
			else:
				if inv_pat:
					return [ gt[:1] + '|' + inverse(gt[2]) + gt[3:] for gt in gts ]
				else:
					return [ gt[:1] + '|' + gt[2:] for gt in gts ]
		
		def inverse_parents_gts(record, inv_mat, inv_pat):
			# both must be hetero
			if inv_mat:
				gt = record.v[9]
				record.v[9] = gt[2] + '|' + gt[0] + gt[3:]
			if inv_pat:
				gt = record.v[10]
				record.v[10] = gt[2] + '|' + gt[0] + gt[3:]
		
		def modify_gts(gts: list[str], record: Optional[VCFFillableRecord]):
			if record is None:
				return
			
			orig_gts = record.v[11:]
			for inv_mat, inv_pat in product((False, True), repeat=2):
				new_gts = inverse_gts(gts, inv_mat, inv_pat)
				if is_near(new_gts, orig_gts):
					inverse_parents_gts(record, inv_mat, inv_pat)
					for c in range(11, len(record.v)):
						record.v[c] = new_gts[c-11]
					break
			else:
				new_gts = gts
			for c in range(11, len(record.v)):
				record.v[c] = new_gts[c-11]
		
		# どちらから来たか決める
		def select_from(froms: list[int], record1: Optional[VCFRecord],
										  record2: Optional[VCFRecord]) -> int:
			if len(froms) == 1:
				return froms[0]
			
			# どちらかは0でない前提
			if froms[0] == 0:
				return froms[1]
			elif froms[1] == 0:
				return froms[0]
			else:	# 両側にRecordがある
				# 近い方を選ぶ
				if record is None or record1 is None or record2 is None:
					# ここには来ないはず
					return 0
				elif record.pos() * 2 < record1.pos() + record2.pos():
					return froms[0]
				else:
					return froms[1]
		
		def modify_parents_type(record):
			if record is not None:
				record.modify_parents_type()
		
		record, mat_record1, mat_record2, pat_record1, pat_record2 = \
															record_set.records()
		if record is None:
			return
		new_gts = record.v[11:]
		for c, gt, mat_gt1, mat_gt2, pat_gt1, pat_gt2 in zip(count(11),
											*map(gts, record_set.records())):
			prev_mat_from = record_set.from_which_chrom_prev_mat(mat_gt1)
			next_mat_from = record_set.from_which_chrom_next_mat(mat_gt2)
			prev_pat_from = record_set.from_which_chrom_prev_pat(pat_gt1)
			next_pat_from = record_set.from_which_chrom_next_pat(pat_gt2)
			if ((prev_mat_from == 0 and next_mat_from == 0) or
					(prev_pat_from == 0 and next_pat_from == 0)):
				new_gts[c-11] = record.v[c]
				continue
			
			mat_froms = unique_list(prev_mat_from, next_mat_from)
			pat_froms = unique_list(prev_pat_from, next_pat_from)
			pairs = [ x for x in product(mat_froms, pat_froms)
										if x[0] != 0 and x[1] != 0 ]
			mat_from, pat_from = select_pair(pairs, gt)
			res = record.v[c][3:]
			if (mat_from, pat_from) == (0, 0):
				if any(mat_from != 0 for mat_from in mat_froms):
					mat_from_nz = select_from(mat_froms,
												mat_record1, mat_record2)
					new_gts[c-11] = record.gt_from_mat(mat_from_nz, c) + res
				elif any(pat_from != 0 for pat_from in pat_froms):
					pat_from_nz = select_from(pat_froms,
												pat_record1, pat_record2)
					new_gts[c-11] = record.gt_from_pat(pat_from_nz, c) + res
				else:
					continue	# phasingしない
			else:
				new_gts[c-11] = record.gt_from_parent(mat_from, pat_from) + res
		
		modify_gts(new_gts, record)		# 両親が間違っていた時の処理
		modify_parents_type(record)
	
	def find_prev_same_type_record(self, i: int, c: int
											) -> Optional[VCFFillableRecord]:
		type = self.records[i].type
		chromosome = self.records[i].v[0]
		for j in range(i - 1, -1, -1):
			record = self.records[j]
			if record.v[0] != chromosome:
				return None
			if record.type == type and record.v[c][:3] != './.':
				return self.records[j]
		else:
			return None
	
	def find_next_same_type_record(self, i: int, c: int
											) -> Optional[VCFFillableRecord]:
		type = self.records[i].type
		chromosome = self.records[i].v[0]
		for j in range(i + 1, len(self.records)):
			record = self.records[j]
			if record.v[0] != chromosome:
				return None
			if record.type == type and record.v[c][:3] != './.':
				return self.records[j]
		else:
			return None
	
	# 家系ごとで./.にしたGenotypeを補完
	def __impute_NA_mat(self, i: int):
		record = self.records[i]
		for c in range(11, len(record.v)):
			if record.v[c][:3] == '.|.':
				self.__impute_NA_mat_each(i, c)
	
	def __impute_NA_mat_each(self, i: int, c: int):
		def select_mat(pairs: list[tuple[int, int]]) -> int:
			if len(pairs) == 1:
				return pairs[0][0]
#			elif prev_record is None or next_record is None:
			elif prev_record is None and next_record is None:
				return 0
			elif next_record is None:
				return prev_mat_from
			elif prev_record is None:
				return next_mat_from
			elif record.pos() * 2 < prev_record.pos() + next_record.pos():
				return prev_mat_from
			else:
				return next_mat_from
		
		record = self.records[i]
		prev_record = self.find_prev_same_type_record(i, c)
		next_record = self.find_next_same_type_record(i, c)
		recordset = RecordSet(record, prev_record, next_record, None, None)
		prev_gt = prev_record.v[c] if prev_record is not None else ''
		next_gt = next_record.v[c] if next_record is not None else ''
		prev_mat_from = recordset.from_which_chrom_prev_mat(prev_gt)
		next_mat_from = recordset.from_which_chrom_next_mat(next_gt)
		mat_froms = unique_list(prev_mat_from, next_mat_from)
		pairs = [ (x, 1) for x in mat_froms if x != 0 ]
		if not pairs:	# 両側ともにmatが無い（まず無い）
			return
		
		mat_from = select_mat(pairs)
		if mat_from == 0:
			return
		
		record.v[c] = record.gt_from_parent(mat_from, 1) + record.v[c][3:]
	
	# 家系ごとで./.にしたGenotypeを補完
	def __impute_NA_pat(self, i: int):
		record = self.records[i]
		for c in range(11, len(record.v)):
			if record.v[c][:3] == './.':
				self.__impute_NA_pat_each(i, c)
	
	def __impute_NA_pat_each(self, i: int, c: int):
		def select_pat(pairs: list[tuple[int, int]]) -> int:
			if len(pairs) == 1:
				return pairs[0][1]
			elif prev_record is None and next_record is None:
				return 0
			elif next_record is None:
				return prev_pat_from
			elif prev_record is None:
				return next_pat_from
			elif record.pos() * 2 < prev_record.pos() + next_record.pos():
				return prev_pat_from
			else:
				return next_pat_from
		
		record = self.records[i]
		prev_record = self.find_prev_same_type_record(i, c)
		next_record = self.find_next_same_type_record(i, c)
		recordset = RecordSet(record, None, None, prev_record, next_record)
		prev_gt = prev_record.v[c] if prev_record is not None else ''
		next_gt = next_record.v[c] if next_record is not None else ''
		prev_pat_from = recordset.from_which_chrom_prev_pat(prev_gt)
		next_pat_from = recordset.from_which_chrom_next_pat(next_gt)
		pat_froms = unique_list(prev_pat_from, next_pat_from)
		pairs = [ (1, x) for x in pat_froms if x != 0 ]
		if not pairs:	# 両側ともにpatが無い（まず無い）
			return
		
		pat_from = select_pat(pairs)
		record.v[c] = record.gt_from_parent(1, pat_from) + record.v[c][3:]
	
	def __impute_others(self, i: int):
		def correct(GT: str, c: int, record: VCFFillableRecord) -> str:
			int_gt = record.get_int_gt(c-9)
			if int_gt == -1:
				if GT[0] != '.':
					return GT[:2] + str(record.pos() % 2)
				elif GT[2] != '.':
					return str(record.pos() % 2) + GT[1:]
				else:
					return (str(record.pos() % 2) + '|' +
									str((record.pos() >> 1) % 2))
			elif int_gt == 0:
				if GT[0] == '1':
					return '1|0'
				elif GT[2] == '1':
					return '0|1'
				else:
					return '0|0'
			elif int_gt == 1:
				if GT[0] == '1':
					return '1|0'
				elif GT[2] == '1':
					return '0|1'
				else:
					# どちらにすべきか判別がつかないので適当に選ぶ
					return '0|1' if int(record.v[1]) % 2 == 0 else '1|0'
			else:
				if GT[0] == '0':
					return '0|1'
				elif GT[2] == '0':
					return '1|0'
				else:
					return '1|1'
		
		record = self.records[i]
		cs = [ c for c in range(11, len(record.v))
				if record.v[c][1] == '/' or record.get_int_gt(c-9) == -1 ]
		if not cs:
			return
		
		mat_homo = record.is_homo(0)
		pat_homo = record.is_homo(1)
		for c in cs:
			mat_from = 1 if mat_homo else self.__find_mat_from(i, c)
			pat_from = 1 if pat_homo else self.__find_pat_from(i, c)
			GT = record.gt_from_parent(mat_from, pat_from)
			if '.' in GT:
				GT = correct(GT, c, record)
			record.set_GT(c-9, GT)
	
	def __find_mat_from(self, i: int, c: int) -> int:
		return self.__select_from(self.__find_prev_mat_from(i, c),
									self.__find_next_mat_from(i, c), i)
	
	def __find_pat_from(self, i: int, c: int) -> int:
		return self.__select_from(self.__find_prev_pat_from(i, c),
									self.__find_next_pat_from(i, c), i)
	
	def __select_from(self, f1: tuple[int, int],
							f2: tuple[int, int], i: int) -> int:
		i1, from1 = f1
		i2, from2 = f2
		if from1 == 0 and from2 == 0:
			# 前後がないとき乱数的に決める
			r0 = self.records[i]
			return r0.pos() % 2 + 1
		
		if from1 == from2:
			return from1
		elif from2 == 0:
			return from1
		elif from1 == 0:
			return from2
		else:
			# 最後は物理距離で決める
			r0 = self.records[i]
			r1 = self.records[i1]
			r2 = self.records[i2]
			if r0.pos() * 2 <= r2.pos() + r1.pos():
				return from1
			else:
				return from2
	
	def __find_prev_mat_from(self, i: int, c: int) -> tuple[int, int]:
		for k in range(i - 1, -1, -1):
			from1 = self.records[k].mat_from(c)
			if from1 != 0:
				return (k, from1)
		return (-1, 0)
	
	def __find_next_mat_from(self, i: int, c: int) -> tuple[int, int]:
		for k in range(i + 1, len(self.records)):
			from1 = self.records[k].mat_from(c)
			if from1 != 0:
				return (k, from1)
		return (-1, 0)
	
	def __find_prev_pat_from(self, i: int, c: int) -> tuple[int, int]:
		for k in range(i - 1, -1, -1):
			from1 = self.records[k].pat_from(c)
			if from1 != 0:
				return (k, from1)
		return (-1, 0)
	
	def __find_next_pat_from(self, i: int, c: int) -> tuple[int, int]:
		for k in range(i + 1, len(self.records)):
			from1 = self.records[k].pat_from(c)
			if from1 != 0:
				return (k, from1)
		return (-1, 0)
	
	def remove_parents(self) -> VCFSmall:
		samples: list[str] = self.samples[2:]
		records: list[VCFRecord] = []
		for record in self.records:
			v = record.v[:9] + record.v[11:]
			r = VCFRecord(v, samples)
			records.append(r)
		header = self.trim_header(samples)
		return VCFSmall(header, records)
	
	@staticmethod
	def fill(vcfs: list[VCFHeteroHomo],
					records: list[VCFImpFamilyRecord]) -> VCFFillable:
		merged_records = VCFFillable.merge_records(vcfs, records)
		vcf: VCFFillable = VCFFillable(vcfs[0].header, merged_records)
		vcf.modify()
		return vcf
	
	@staticmethod
	def merge_records(vcfs: list[VCFHeteroHomo],
						records: list[VCFImpFamilyRecord],
						) -> list[VCFFillableRecord]:
		all_records = [ VCFFillableRecord.convert(record) for record in
						chain((r for vcf in vcfs for r in vcf.records), records)
		]
		all_records.sort(key=lambda record: record.index)
		return all_records
	
	@staticmethod
	def collect_records(vcfs: list[VCFFillable]
									) -> list[list[VCFFillableRecord]]:
		rss: list[list[VCFFillableRecord]] = []
		if not vcfs:
			return rss
		
		for i in range(len(vcfs[0])):
			rss.append([ vcf.records[i] for vcf in vcfs ])
		return rss
	
	@staticmethod
	def integrate_samples(sss: list[list[str]], orig_samples: list[str]
							) -> tuple[list[str], list[list[tuple[int, int]]]]:
		dic = defaultdict(list)
		for i, samples in enumerate(sss):
			for j, sample in enumerate(samples):
				dic[sample].append((i, j))
		samples = [ sample for sample in orig_samples if sample in dic ]
		pos_samples = [ dic[sample] for sample in samples ]
		return (samples, pos_samples)
	
	# 重複したサンプルが一つになるようにVCFを統合する
	@staticmethod
	def integrate(vcf: VCFFillable, rss: list[list[VCFFillableRecord]],
										orig_samples: list[str]) -> VCFSmall:
		sss = [ r.samples for r in rss[0] ]
		samples, pos_samples = VCFFillable.integrate_samples(sss, orig_samples)
		records: list[VCFRecord] = []
		for rs in rss:
			records.append(VCFFillableRecord.integrate(rs, samples, pos_samples))
		header = vcf.trim_header(records[0].samples)
		new_vcf = VCFSmall(header, [])
		new_vcf.records.extend(records)
		return new_vcf
	
	@staticmethod
	def merge(vcfs: list[VCFFillable], orig_samples: list[str]) -> VCFSmall:
		rss = VCFFillable.collect_records(vcfs)
		return VCFFillable.integrate(vcfs[0], rss, orig_samples)


__all__ = ['VCFFillableRecord', 'VCFFillable']
