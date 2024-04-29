#include "../include/LargeFamily.h"
#include "../include/VCFHeteroHomo.h"
#include "../include/VCFFillable.h"
#include "../include/VCFImpFamily.h"
#include "../include/VCFJunkRecord.h"
#include "../include/VCFHomoHomo.h"
#include "../include/VCFHeteroHeteroLite.h"
#include "../include/KnownFamily.h"
#include "../include/Map.h"
#include "../include/option.h"
#include "../include/common.h"

using namespace std;


//////////////////// LargeFamily ////////////////////

VCFHeteroHomoRecord *LargeFamily::create_heterohomo_record(const STRVEC& v,
										const KnownFamily *family, size_t i,
										WrongType wrong_type, ParentComb pc) {
	const STRVEC&	samples = family->get_samples();
	const int	total_int_gt = pc == ParentComb::P00x01 ? 1 : 3;
	if(!family->is_mat_known()) {
		const int	pat_int_gt = Genotype::get_int_gt(v[10]);
		const int	mat_int_gt = total_int_gt - pat_int_gt;
		if(mat_int_gt < 0 || 2 < mat_int_gt)
			wrong_type = WrongType::MODIFIABLE;
		else {
			STRVEC	w(v.size());
			std::copy(v.begin(), v.end(), w.begin());
			w[9] = Genotype::int_to_gt(mat_int_gt);
			return new VCFHeteroHomoRecord(w, samples, i, wrong_type, pc);
		}
	}
	else if(!family->is_pat_known()) {
		const int	mat_int_gt = Genotype::get_int_gt(v[9]);
		const int	pat_int_gt = total_int_gt - mat_int_gt;
		if(pat_int_gt < 0 || 2 < pat_int_gt)
			wrong_type = WrongType::MODIFIABLE;
		else {
			STRVEC	w(v.size());
			std::copy(v.begin(), v.end(), w.begin());
			w[10] = Genotype::int_to_gt(pat_int_gt);
			return new VCFHeteroHomoRecord(w, samples, i, wrong_type, pc);
		}
	}
	
	return new VCFHeteroHomoRecord(v, samples, i, wrong_type, pc);
}

void LargeFamily::classify_record(size_t i, const VCFFamily *vcf,
								  const TypeDeterminer *td,
								  const KnownFamily *family,
								  vector<VCFHeteroHomoRecord *>& heho_records,
								  vector<VCFImpFamilyRecord *>& other_records) {
	ClassifyRecord	*CR = ClassifyRecord::get_instance();
	VCFFamilyRecord	*record = vcf->get_family_record(i);
	const auto	pair1 = CR->classify(record, td, family->is_one_unknown());
	const ParentComb	pc = pair1.first;
	const WrongType	wrong_type = pair1.second;
	const STRVEC&	v = record->get_v();
	const STRVEC&	samples = record->get_samples();
	if(pc == ParentComb::PNA) {		// no candidates
		other_records[i] = new VCFJunkRecord(v, samples, i, wrong_type);
	}
	else if(TypeDeterminer::is_homohomo(pc)) {
		auto	*r_ = new VCFHomoHomoRecord(v, samples, i, wrong_type, pc);
		auto	*r = r_->impute();
		delete r_;
		other_records[i] = r;
	}
	else if(TypeDeterminer::is_heterohomo(pc)) {
		heho_records[i] = create_heterohomo_record(v, family,
													i, wrong_type, pc);
	}
	else {
		other_records[i] = new VCFHeteroHeteroLiteRecord(v, samples, i,
															wrong_type, pc);
	}
}

void LargeFamily::classify_records_in_thread(void *config) {
	auto	*c = (ConfigThreadClassify *)config;
	const size_t	n = c->vcf->size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		classify_record(i, c->vcf, c->td, c->family,
							c->heho_records, c->other_records);
	}
}

pair<vector<VCFHeteroHomoRecord *>, vector<VCFImpFamilyRecord *>>
LargeFamily::classify_records(const VCFFamily *vcf,
							  const KnownFamily *family, const Option *option) {
	const size_t	N = vcf->size();
	vector<VCFHeteroHomoRecord *>	heho_records(N, NULL);
	vector<VCFImpFamilyRecord *>	other_records(N, NULL);
	ClassifyRecord	*CR = ClassifyRecord::get_instance();
	const auto	*td = CR->get_TypeDeterminer(vcf->num_samples()-2,
														option->ratio);
	const int	T = vcf->size() < 100 ? 1 : option->num_threads;
	
	vector<ConfigThreadClassify *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThreadClassify(i, T, vcf, td, family,
											  heho_records, other_records);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
			(void *(*)(void *))&classify_records_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		classify_records_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	
	return std::make_pair(heho_records, other_records);
}

void LargeFamily::create_vcf_in_thread(void *config) {
	auto	*c = (ConfigThreadCreate *)config;
	const VCFSmall	*orig_vcf = c->orig_vcf;
	for(size_t i = c->first; i < c->families.size(); i += c->num_threads) {
		const KnownFamily	*family = c->families[i];
		c->results[i] = VCFFamily::create(orig_vcf, family->get_samples());
	}
}

vector<VCFFamily *> LargeFamily::create_family_vcfs(
								const VCFSmall *orig_vcf,
								const vector<const KnownFamily *>& families,
								int num_threads) {
	const size_t	N = families.size();
	vector<VCFFamily *>	vcfs_family(N, NULL);
	const int	T = std::min((int)N, num_threads);
	vector<ConfigThreadCreate *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThreadCreate(i, T, orig_vcf,
											families, vcfs_family);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
			(void *(*)(void *))&create_vcf_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		create_vcf_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	
	return vcfs_family;
}

vector<pair<string, vector<pair<size_t, size_t>>>>
LargeFamily::collect_same_parent_families(
								const vector<const KnownFamily *>& families) {
	// { parent: [(family index, parent index)] }
	map<string, vector<pair<size_t, size_t>>>	dic;
	for(size_t i = 0; i < families.size(); ++i) {
		dic[families[i]->get_mat()].push_back(make_pair(i, 0));
		dic[families[i]->get_pat()].push_back(make_pair(i, 1));
	}
	
	vector<pair<string, vector<pair<size_t, size_t>>>>	v;
	for(auto p = dic.begin(); p != dic.end(); ++p) {
		if(p->second.size() > 1)
			v.push_back(*p);
	}
	return v;
}

vector<vector<VCFImpFamilyRecord *>> LargeFamily::collect_same_position_records(
				const vector<tuple<vector<VCFHeteroHomoRecord *>,
								   vector<VCFImpFamilyRecord *>, int>>& rs) {
	const size_t	vcf_size = get<0>(rs[0]).size();
	vector<vector<VCFImpFamilyRecord *>>	recordss(vcf_size);
	for(auto p = rs.begin(); p != rs.end(); ++p) {
		const auto&	rs1 = get<0>(*p);
		for(auto q = rs1.begin(); q != rs1.end(); ++q) {
			auto	*r1 = *q;
			if(r1 != NULL && r1->is_right())
				recordss[r1->get_index()].push_back(r1);
		}
		const auto&	rs2 = get<1>(*p);
		for(auto q = rs2.begin(); q != rs2.end(); ++q) {
			auto	*r2 = *q;
			if(r2 != NULL && r2->is_homohomo())
				recordss[r2->get_index()].push_back(r2);
		}
	}
	return recordss;
}

void LargeFamily::modify_00x11_each(
				const vector<tuple<vector<VCFHeteroHomoRecord *>,
								   vector<VCFImpFamilyRecord *>, int>>& rs) {
	auto	records = collect_same_position_records(rs);
	
	for(auto p = records.begin(); p != records.end(); ++p) {
		if(p->size() >= 2)
			VCFImpFamilyRecord::modify_00x11(*p);
	}
}

// When a parent of record with 0/0 x 1/1 parents is one of another family,
// correct if the genotype is different from other Homo x Homo or Hetero x Homo
void LargeFamily::modify_00x11(
					vector<vector<VCFHeteroHomoRecord *>>& heho_recordss,
				 	vector<vector<VCFImpFamilyRecord *>>& other_recordss,
				 	const vector<const KnownFamily *>& families) {
	auto	fams = collect_same_parent_families(families);
	for(auto p = fams.begin(); p != fams.end(); ++p) {
		auto	v = p->second;
		vector<tuple<vector<VCFHeteroHomoRecord *>,
					 vector<VCFImpFamilyRecord *>, int>>	rs;
		for(auto q = v.begin(); q != v.end(); ++q) {
			const size_t	fam_index = q->first;
			const size_t	p_index = q->second;
			rs.push_back(make_tuple(heho_recordss[fam_index],
									other_recordss[fam_index], p_index));
		}
		modify_00x11_each(rs);
	}
}

// divide records into hetero x homo and other types
pair<vector<vector<VCFHeteroHomoRecord *>>,
	 vector<vector<VCFImpFamilyRecord *>>>
LargeFamily::divide_vcf_into_record_types(
								const vector<VCFFamily *>& family_vcfs,
								const vector<const KnownFamily *>& families,
								const Option *option) {
	vector<vector<VCFHeteroHomoRecord *>>	heho_recordss;
	vector<vector<VCFImpFamilyRecord *>>	other_recordss;
	for(size_t i = 0; i < family_vcfs.size(); ++i) {
		auto	*vcf = family_vcfs[i];
		auto	*family = families[i];
		auto	pair1 = classify_records(vcf, family, option);
		auto	heho_records = pair1.first;
		auto	other_records = pair1.second;
		heho_recordss.push_back(heho_records);
		other_recordss.push_back(other_records);
	}
	return make_pair(heho_recordss, other_recordss);
}

void LargeFamily::fill_in_thread(void *config) {
	const auto	*c = (ConfigThreadFill *)config;
	const size_t	n = c->vcfss_heho.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		const auto&	vcfs = c->vcfss_heho[i];
		const auto&	records = c->other_recordss[i];
		c->results[i] = VCFFillable::fill(vcfs, records,
											c->num_threads_in_family);
	}
}

vector<VCFFillable *> LargeFamily::fill_vcf(
					const map<string, vector<VCFHeteroHomo *>>& dic_vcfs,
					const vector<vector<VCFImpFamilyRecord *>>& other_recordss,
					const vector<const KnownFamily *>& families,
					int num_threads) {
	// collect vcf by family index
	const size_t	N = families.size();
	map<pair<string, string>, size_t>	family_indices;
	for(size_t i = 0; i < N; ++i)
		family_indices.insert(make_pair(families[i]->parents(), i));
	vector<vector<VCFHeteroHomo *>>	vcfss_heho(N);
	for(auto p = dic_vcfs.begin(); p != dic_vcfs.end(); ++p) {
		const auto&	vcfs = p->second;
		for(auto q = vcfs.begin(); q != vcfs.end(); ++q) {
			auto	*vcf = *q;
			const size_t	index = family_indices[vcf->parents()];
			vcfss_heho[index].push_back(vcf);
		}
	}
	
	// insert and correct other records
	vector<VCFFillable *>	vcfs_filled(N);
	const int	T = std::min(num_threads, (int)N);
	const int	T_in_family = std::max(2, (T+(int)N-1)/(int)N);
	vector<ConfigThreadFill *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThreadFill(i, T, vcfss_heho, other_recordss,
													T_in_family, vcfs_filled);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
					(void *(*)(void *))&fill_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		fill_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	return vcfs_filled;
}

void LargeFamily::compress_records(vector<VCFImpFamilyRecord *>& others) {
	auto	p = others.begin();
	for(auto q = others.begin(); q != others.end(); ++q) {
		if(*q != NULL) {
			*p = *q;
			++p;
		}
	}
	others.erase(p, others.end());
}

void LargeFamily::clean_in_thread(void *config) {
	auto	*c = (ConfigThreadClean *)config;
	const int	n = (int)c->families.size();
	const int	T = c->num_threads;
	const int	T_in_family = std::max(2, (T+n-1)/n);
	for(int i = c->first; i < n; i += T) {
		const auto&	records = c->recordss[i];
		const auto&	samples = c->families[i]->get_samples();
		auto	header = c->orig_vcf->trim_header(samples);
		auto	p = VCFHeteroHomo::clean_vcfs(records, header, samples,
													c->geno_map, T_in_family);
		c->vcfss[i] = p.first;
		c->unused_recordss[i] = p.second;
	}
}

VCFSmall *LargeFamily::correct_large_family_VCFs(
									const VCFSmall *orig_vcf,
									const vector<const KnownFamily *>& families,
									const Map& geno_map, const Option *option) {
	const int	num_threads = option->num_threads;
	// create a VCF for each large family
	const auto	family_vcfs = create_family_vcfs(orig_vcf,
												families, num_threads);
	// classify records
	auto	p = divide_vcf_into_record_types(family_vcfs, families, option);
	auto&	heho_recordss = p.first;
	auto&	other_recordss = p.second;
	
	// We have to deal with it all together
	// because we coordinate between families.
	modify_00x11(heho_recordss, other_recordss, families);
	
	// impute VCFs by family
	const size_t	LN = families.size();
	vector<vector<VCFHeteroHomo *>>	vcfss(LN);
	vector<vector<VCFHeteroHomoRecord *>>	unused_recordss(LN);
	const int	T = min((int)LN, num_threads);
	vector<ConfigThreadClean *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThreadClean(i, T, heho_recordss,
										   orig_vcf, families,
										   geno_map, vcfss, unused_recordss);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
					(void *(*)(void *))&clean_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		clean_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	
	for(size_t i = 0; i < families.size(); ++i) {
		auto&	others = other_recordss[i];
		const auto&	unused = unused_recordss[i];
		compress_records(others);
		others.insert(others.end(), unused.begin(), unused.end());
	}
	
	Common::delete_all(family_vcfs);
	
	// collect VCFs by parent and decide phasing
	map<string, vector<VCFHeteroHomo *>>	dic_vcfs;
	for(size_t i = 0; i < families.size(); ++i) {
		auto	vcfs = vcfss[i];
		const Family	*family = families[i];
		for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
			auto	*vcf = *p;
			if(vcf->size() == 0) {
				// Add a VCF to a dic_vcfs even if it is empty,
				// because if no VCF exists, you won't know the family
				dic_vcfs[family->get_mat()].push_back(vcf);
				dic_vcfs[family->get_pat()].push_back(vcf);
			}
			else if(vcf->is_mat_hetero())
				dic_vcfs[family->get_mat()].push_back(vcf);
			else
				dic_vcfs[family->get_pat()].push_back(vcf);
		}
	}
	
	for(auto p = dic_vcfs.begin(); p != dic_vcfs.end(); ++p) {
		auto	vcfs = p->second;
		VCFHeteroHomo::inverse_phases(vcfs);
	}
	
	auto	vcfs_filled = fill_vcf(dic_vcfs, other_recordss,
										families, option->num_threads);
	for(size_t i = 0; i < families.size(); ++i) {
		Common::delete_all(vcfss[i]);
		Common::delete_all(other_recordss[i]);
	}
	auto	vcf1 = VCFFillable::merge(vcfs_filled, orig_vcf->get_samples());
	for(auto p = vcfs_filled.begin(); p != vcfs_filled.end(); ++p) {
		delete *p;
	}
	return vcf1;
}
