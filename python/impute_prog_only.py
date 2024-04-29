# coding: utf-8
# impute_prog_only.py

from __future__ import annotations
from itertools import *
import sys

from VCF import *
from VCFFamily import VCFFamily
from VCFImpFamily import FillType, VCFImpFamilyRecord
from VCFHeteroHomo import VCFHeteroHomoRecord
from VCFHeteroHomoPP import *
from VCFHomoHomo import VCFHomoHomoRecord
from VCFHeteroHeteroLite import VCFHeteroHeteroLiteRecord
from VCFJunkRecord import VCFJunkRecord
from VCFFillableRecord import VCFFillableRecord
from VCFFillable import VCFFillable
import SmallFamily
from pedigree import PedigreeTable, Family
from SampleManager import SampleManager
import ClassifyRecord as CR
from TypeDeterminer import ParentComb, TypeDeterminer
from Map import *
from option import *
from common import *


#################### process ####################

def classify_record(i: int, record: VCFRecord, td: TypeDeterminer,
							heho_records: list[Optional[VCFHeteroHomoRecord]],
							other_records: list[Optional[VCFImpFamilyRecord]]):
	samples = record.samples
	wrong_type, pair = CR.classify_record(record, td, False)
	v = record.v
	if pair.is_homohomo():
		record_ = VCFHomoHomoRecord(v, samples, i, wrong_type, pair)
		other_records[i] = record_.impute()
	elif pair.is_heterohomo():
		heho_records[i] = VCFHeteroHomoRecord(v, samples, i,
													wrong_type, pair)
	elif pair == ParentComb.P01x01:
		other_records[i] = VCFHeteroHeteroLiteRecord(v, samples, i, wrong_type)
	else:		# 候補が無い
		other_records[i] = VCFJunkRecord(v, samples, i, wrong_type)

def classify_records(vcf: VCFFamily, option: Option
									) -> tuple[list[VCFHeteroHomoRecord],
												 list[VCFImpFamilyRecord]]:
	td = CR.get_typedeterminer(vcf.num_samples()-2, option.ratio)
	records1: list[Optional[VCFHeteroHomoRecord]] = [None] * len(vcf)
	records2: list[Optional[VCFImpFamilyRecord]] = [None] * len(vcf)
	for i, record in enumerate(vcf.records):
		classify_record(i, record, td, records1, records2)
	
	heho_records: list[VCFHeteroHomoRecord] = []
	other_records: list[VCFImpFamilyRecord] = []
	for i in range(len(vcf)):
		r1 = records1[i]
		r2 = records2[i]
		if r1 is not None:
			heho_records.append(r1)
		elif r2 is not None:
			other_records.append(r2)
	
	return (heho_records, other_records)

def merge_vcf(rss: list[list[VCFFillableRecord]],
					header: list[list[str]], gmap: Map) -> VCFHeteroHomoPP:
	rs = rss[0] + rss[1] + rss[2] + rss[3]
	rs.sort(key=lambda r: r.pos())
	return VCFHeteroHomoPP(header, rs, gmap)

def impute_prog_vcf_chr(parent_vcf: VCFSmall, prog_vcf: VCFSmall,
						gmap: Map, option: Option) -> VCFHeteroHomoPP:
	print('chr: %s %d records' % (parent_vcf.records[0].chrom(),
													len(parent_vcf)))
	sys.stdout.flush()
	# parent_vcfは両親のみ、prog_vcfは後代のみという前提
	samples = parent_vcf.samples + prog_vcf.samples
	merged_vcf = VCFHeteroHomoPP.merge(parent_vcf, prog_vcf,
											samples, gmap, option)
	rss = VCFHeteroHomoPP.classify_records(merged_vcf.records)
	header = parent_vcf.trim_header(samples)
	mat_vcf = VCFHeteroHomoPP(header, rss[FillType.MAT.value], gmap)
	pat_vcf = VCFHeteroHomoPP(header, rss[FillType.PAT.value], gmap)
	mat_vcf.impute()
	pat_vcf.impute()
	merged_vcf = merge_vcf(rss, header, gmap)
	merged_vcf.fill()
	return merged_vcf

__all__ = ['impute_prog_vcf_chr']
