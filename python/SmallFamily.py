from __future__ import annotations

# coding: utf-8
# SmallFamily.py

from itertools import *
from collections import defaultdict, Counter
import random
import csv
import sys
from multiprocessing import Pool

from typing import Dict, List, Tuple, Set, Optional

from VCF import *
from VCFFamily import VCFFamily, VCFFamilyBase, VCFFamilyRecord
from pedigree import PedigreeTable, Family
from VCFHeteroHomoPP import *
from VCFOneParentPhased import VCFOneParentPhased
from VCFProgenyPhased import VCFProgenyPhased
from VCFIsolated import VCFIsolated
from SampleManager import *
from Map import *
from option import *
from common import *


def impute_vcf_by_parents_core(orig_vcf: VCFSmall, merged_vcf: VCFSmallBase,
								families: list[Family], gmap: Map) -> VCFSmall:
	vcfs: list[VCFSmallBase] = []
	for family in families:
		family_vcf = VCFHeteroHomoPP.impute_by_parents(
								orig_vcf, merged_vcf, family.samples(), gmap)
		vcf = family_vcf.remove_parents()
		vcfs.append(vcf)
	
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	return new_vcf

def impute_vcf_by_parents(orig_vcf: VCFSmall, merged_vcf: VCFSmallBase,
							sample_man: SampleManager,
							gmap: Map, num_threads: int) -> Optional[VCFSmall]:
	families = sample_man.extract_small_families()
	if not families:
		return None
	
	new_imputed_vcf = impute_vcf_by_parents_core(orig_vcf, merged_vcf,
															families, gmap)
	merged_vcf = VCFSmall.join([merged_vcf, new_imputed_vcf], orig_vcf.samples)
	sample_man.set(new_imputed_vcf.samples)
	return merged_vcf

def impute_vcf_by_parent_core(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
										families: list[Family], gmap: Map,
										sample_man: SampleManager,
										num_threads: int) -> VCFSmall:
	# collect not phased parents
	samples: list[str] = [ parent for family in families
								  for parent in family.parents()
								  if not sample_man.is_imputed(parent) ]
	
	# phase not phased parents
	parents_vcf = impute_isolated_samples(orig_vcf, merged_vcf, sample_man,
													samples, gmap, num_threads)
	
	# merge vcfs
	new_merged_vcf = VCFSmall.join([merged_vcf, parents_vcf], orig_vcf.samples)
	
	# impute progenies
	new_vcf = impute_vcf_by_parents_core(orig_vcf, new_merged_vcf,
														families, gmap)
	if new_vcf is None:
		return new_merged_vcf	# not comes here
	
	# join
	return VCFSmall.join([parents_vcf, new_vcf], orig_vcf.samples)

def impute_vcf_by_parent(orig_vcf: VCFSmall, merged_vcf: VCFSmall, gmap: Map,
										sample_man: SampleManager,
										num_threads: int) -> Optional[VCFSmall]:
	families = sample_man.extract_single_parent_phased_families()
	if not families:
		return None
	
	new_imputed_vcf = impute_vcf_by_parent_core(orig_vcf, merged_vcf, families,
												gmap, sample_man, num_threads)
	merged_vcf = VCFSmall.join([merged_vcf, new_imputed_vcf],
												orig_vcf.samples)
	sample_man.set(new_imputed_vcf.samples)
	return merged_vcf

def impute_one_parent_vcf_core(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
										families: list[Family], gmap: Map,
										sample_man: SampleManager,
										num_threads: int) -> VCFSmall:
	references = sample_man.collet_references()
	ref_vcf = merged_vcf.extract_samples(references)
	
	vcfs: list[VCFSmallBase] = []
	for family in families:
		is_mat_phased = sample_man.is_imputed(family.mat)
		opp_vcf = VCFOneParentPhased.create(family.samples(), is_mat_phased,
											merged_vcf, orig_vcf, gmap, ref_vcf)
		opp_vcf.impute()
		vcfs.append(opp_vcf)
	
	samples = []	# 今回phasingしたsampleを集める
	for family in families:
		ss = [ s for s in family.samples()
					if not sample_man.is_imputed(s) and s != '0' ]
		samples.extend(ss)
		sample_man.set(ss)
	
	return VCFSmall.join(vcfs, samples)

def impute_one_parent_vcf(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
										gmap: Map, sample_man: SampleManager,
										num_threads: int) -> Optional[VCFSmall]:
	families = sample_man.extract_phased_and_unknown_parents_family()
	if not families:
		return None
	
	new_imputed_vcf = impute_one_parent_vcf_core(orig_vcf, merged_vcf, families,
												 gmap, sample_man, num_threads)
	merged_vcf = VCFSmall.join([merged_vcf, new_imputed_vcf],
												orig_vcf.samples)
	sample_man.set(new_imputed_vcf.samples)
	return merged_vcf

def impute_vcf_by_progenies_core(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
											families: list[Family], gmap: Map,
											sample_man: SampleManager,
											num_threads: int) -> VCFSmall:
	references = sample_man.collet_references()
	ref_vcf = merged_vcf.extract_samples(references)
	vcfs: list[VCFSmallBase] = []
	for family in families:
		ppi = next(i for i, sample in enumerate(family.samples())
								if sample_man.is_imputed(sample))
		pp_vcf = VCFProgenyPhased.create(orig_vcf, merged_vcf,
									family.samples(), ppi, gmap, ref_vcf)
		pp_vcf.impute()		# impute only parents
		if pp_vcf.num_progenies() == 1:
			vcf4: VCFSmallBase = pp_vcf	# all imputed
		elif family.mat != '0' and family.pat != '0':
			vcf2 = impute_vcf_by_parents(orig_vcf, pp_vcf,
											sample_man, gmap, num_threads)
			if vcf2 is None:
				return merged_vcf	# not comes here
			vcf4 = VCFSmall.join([pp_vcf, vcf2], family.samples())
		else:
			is_mat_phased = family.mat != '0'
			vcf3 = VCFOneParentPhased.create(family.samples(), is_mat_phased,
												pp_vcf, orig_vcf, gmap, ref_vcf)
			vcf3.impute()
			vcf4 = VCFSmall.join([pp_vcf, vcf3], family.samples())
		vcfs.append(vcf4)
	
	samples = []	# collect phased samples
	for family in families:
		ss = [ s for s in family.samples()
							if not sample_man.is_imputed(s) ]
		samples.extend(ss)
		sample_man.set(ss)
	
	return VCFSmall.join(vcfs, samples)

def impute_vcf_by_progenies(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
									gmap: Map, sample_man: SampleManager,
									num_threads: int) -> Optional[VCFSmall]:
	families = sample_man.extract_progenies_phased_families()
	if not families:
		return None
	
	new_imputed_vcf = impute_vcf_by_progenies_core(orig_vcf, merged_vcf,
													families, gmap,
													sample_man, num_threads)
	merged_vcf = VCFSmall.join([merged_vcf, new_imputed_vcf],
												orig_vcf.samples)
	sample_man.set(new_imputed_vcf.samples)
	return merged_vcf

def impute_small_family_VCFs(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
								geno_map: Map, sample_man: SampleManager,
								num_threads: int) -> VCFSmall:
	while True:
		# 両親が補完されているが子どもが少ない家系を補完する
		new_merged_vcf1 = impute_vcf_by_parents(orig_vcf, merged_vcf,
											sample_man, geno_map, num_threads)
		if new_merged_vcf1 is not None:
			merged_vcf = new_merged_vcf1
			continue
		
		# 片親が補完されている家系を補完する
		new_merged_vcf2 = impute_vcf_by_parent(orig_vcf, merged_vcf,
											geno_map, sample_man, num_threads)
		if new_merged_vcf2 is not None:
			merged_vcf = new_merged_vcf2
			continue
		
		new_merged_vcf3 = impute_one_parent_vcf(orig_vcf, merged_vcf,
											geno_map, sample_man, num_threads)
		if new_merged_vcf3 is not None:
			merged_vcf = new_merged_vcf3
			continue
		
		# Impute families whose progenies have been imputed
		new_merged_vcf4 = impute_vcf_by_progenies(orig_vcf, merged_vcf,
											geno_map, sample_man, num_threads);
		if new_merged_vcf4 is not None:
			merged_vcf = new_merged_vcf4
			continue
		
		break
	
	return merged_vcf

def impute_isolated_samples(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
									sample_man: SampleManager,
									samples: list[str],
									gmap: Map, num_threads: int) -> VCFSmall:
	references = sample_man.collet_references()
	# あとでマルチプロセス化するためにphasingすべきsample分割する
	vcfs = VCFIsolated.create(orig_vcf, merged_vcf,
								samples, references, gmap, num_threads)
	new_vcfs: list[VCFSmallBase] = []
	for vcf in vcfs:
		vcf.impute()
		new_vcfs.append(vcf.extract_isolated_samples())
	
	new_vcf = VCFSmall.join(new_vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute_small_family_VCFs', 'impute_isolated_samples']
