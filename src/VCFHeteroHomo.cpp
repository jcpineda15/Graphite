#include <stack>
#include <algorithm>
#include <cmath>
#include <random>
#include <cassert>
#include "../include/common.h"
#include "../include/VCFHeteroHomo.h"
#include "../include/VCFOriginal.h"
#include "../include/Map.h"
#include "../include/Pedigree.h"
#include "../include/Imputer.h"
#include "../include/inverse_graph.h"
#include "../include/option.h"

using namespace std;


//////////////////// VCFHeteroHomoRecord ////////////////////

FillType VCFHeteroHomoRecord::get_fill_type() const {
	if(this->is_imputable())
		return this->is_mat_hetero() ? FillType::MAT : FillType::PAT;
	else
		return FillType::IMPUTABLE;
}

vector<int> VCFHeteroHomoRecord::genotypes_from_hetero_parent() const {
	const vector<int>	gts = this->get_int_gts();
	const int	homo_gt = this->is_mat_hetero() ? pat_int_gt() : mat_int_gt();
	vector<int>	which_froms;
	for(auto p = gts.begin() + 2; p != gts.end(); ++p) {
		const int	gt = homo_gt == 0 ? *p : *p - 1;
		which_froms.push_back(gt != 0 && gt != 1 ? -1 : gt);
	}
	return which_froms;
}

void VCFHeteroHomoRecord::set_haplo(int h) {
	const int	hetero_col = this->is_mat_hetero() ? 9 : 10;
	const int	homo_col = hetero_col == 9 ? 10 : 9;
	this->v[hetero_col] = h == 0 ? "0|1" : "1|0";
	this->v[homo_col] = this->v[homo_col].c_str()[0] == '0' ? "0|0" : "1|1";
	
	const auto	gts = this->genotypes_from_hetero_parent();
	for(size_t i = 0; i < gts.size(); ++i)
		this->which_comes_from[i] = gts[i] == -1 ? -1 : (gts[i] == h ? 0 : 1);
}

void VCFHeteroHomoRecord::set_int_gt_by_which_comes_from(int w, int i) {
	this->which_comes_from[i] = w;
	const string&	mat_gt = this->v[9];
	const string&	pat_gt = this->v[10];
	if(this->is_mat_hetero()) {
		const string	homo_gt = pat_gt.substr(0, 1);
		this->set_GT(i+2, mat_gt.substr(w*2, 1) + '|' + homo_gt);
	}
	else {
		const string	homo_gt = mat_gt.substr(0, 1);
		this->set_GT(i+2, homo_gt + '|' + pat_gt.substr(w*2, 1));
	}
}

//////////////////// VCFHeteroHomo ////////////////////

VCFHeteroHomo::VCFHeteroHomo(const vector<STRVEC>& h, const STRVEC& s,
							vector<VCFHeteroHomoRecord *> rs, const Map& m) :
						VCFBase(h, s), VCFSmallBase(),
						VCFFamilyBase(), VCFMeasurable(m), records(rs) { }

VCFHeteroHomo::~VCFHeteroHomo() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

InvGraph VCFHeteroHomo::make_graph(double max_dist) const {
	vector<vector<int>>	gtss;
	for(auto p = this->records.begin(); p != this->records.end(); ++p) {
		gtss.push_back((*p)->genotypes_from_hetero_parent());
	}
	
	vector<double>	cMs;
	const size_t	L = this->size();
	for(size_t i = 0; i < L; ++i) {
		cMs.push_back(this->record_cM(i));
	}
	
	InvGraph	graph;
	// isolated nodes can not be represented without pre-registering the key
	for(size_t k = 0U; k < L; ++k)
		graph[k];
	for(size_t k = 0; k < L; ++k) {
		for(size_t l = k + 1; l < L; ++l) {
			const double	cM = cMs[l] - cMs[k];
			// don't see a connection
			// if the records are more than 10 cM apart
			if(cM > 10.0)
				break;
			const auto	p = distance(gtss[k], gtss[l]);
			const double	dist = p.first;
			const bool		inv = p.second;
			if(dist <= max_dist) {
				graph[k].push_back(make_tuple(l, dist, inv));
				graph[l].push_back(make_tuple(k, dist, inv));
			}
		}
	}
	return graph;
}

// Since beta distribution is used, dist is a double
pair<double, bool> VCFHeteroHomo::distance(const vector<int>& gts1,
											const vector<int>& gts2) {
	const int	N = (int)gts1.size();
	int	counter_right = 0;
	int	counter_diff = 0;
	int	counter_NA = 0;
	for(size_t i = 0; i < gts1.size(); ++i) {
		const int	int_gt1 = gts1[i];
		const int	int_gt2 = gts2[i];
		if(int_gt1 == -1 || int_gt2 == -1)
			counter_NA += 1;
		else if(int_gt1 == int_gt2)
			counter_right += 1;
		else
			counter_diff += 1;
	}
	
	if(counter_right >= counter_diff)
		return make_pair(dist_with_NA(counter_right, counter_NA, N), false);
	else
		return make_pair(dist_with_NA(counter_diff, counter_NA, N), true);
}

// beta distribution
double VCFHeteroHomo::dist_with_NA(int right, int counter_NA, int N) {
	const int	diff = N - right - counter_NA;
	const double	diff_ratio = (diff + 1.0) / (right + diff + 2.0);
	return diff + counter_NA * diff_ratio;
}

// decide heterozygous parent haplotypes
vector<int> VCFHeteroHomo::make_parent_haplotypes(const InvGraph& graph) const {
	const InvGraph	tree = graph.minimum_spanning_tree();
	auto	vs = tree.collect_nodes();
	std::sort(vs.begin(), vs.end());
	vector<int>	haplo = create_haplotype(vs[0], tree);
	
	vector<int>	haplo2;
	for(auto p = haplo.begin(); p != haplo.end(); ++p) {
		if(*p != -1)
				haplo2.push_back(*p);
	}
	return haplo2;
}

vector<int> VCFHeteroHomo::create_haplotype(size_t v0,
											const InvGraph& tree) const {
	const size_t	L = size();
	vector<int>	haplo(L, -1);
	haplo[v0] = 0;
	const auto	edges = tree.walk(v0);
	for(auto p = edges.begin(); p != edges.end(); ++p) {
		const size_t	v1 = get<0>(*p);
		const size_t	v2 = get<1>(*p);
		const bool		inv = get<2>(*p);
		if(inv)
			haplo[v2] = haplo[v1] == 0 ? 1 : 0;		// inverse
		else
			haplo[v2] = haplo[v1];
	}
	return haplo;
}

VCFHeteroHomo *VCFHeteroHomo::make_subvcf(const InvGraph& graph) const {
	const vector<int>	haplo = make_parent_haplotypes(graph);
	vector<int>	indices;
	for(auto p = graph.begin(); p != graph.end(); ++p)
		indices.push_back(p->first);
	std::sort(indices.begin(), indices.end());
	vector<VCFHeteroHomoRecord *>	records;
	for(size_t i = 0; i < haplo.size(); ++i) {
		VCFHeteroHomoRecord	*record = this->records[indices[i]];
		record->set_haplo(haplo[i]);
		records.push_back(record);
	}
	return new VCFHeteroHomo(this->header, this->samples,
										records, this->get_map());
}

bool operator <(const InvGraph& g1, const InvGraph& g2) {
	return g1.size() < g2.size();
}

pair<vector<VCFHeteroHomo *>, vector<VCFHeteroHomoRecord *>>
		VCFHeteroHomo::determine_haplotype(const OptionImpute *option) const {
	const double	max_dist = std::min((double)option->max_dist,
										this->num_progenies() * 0.1);
	const InvGraph	graph = this->make_graph(max_dist);
	const vector<InvGraph>	subgraphs = graph.divide_into_connected();
	
	// keep only large graphs and collect records of small graphs
	vector<InvGraph>	large_graphs;
	vector<InvGraph>	small_graphs;
	for(auto p = subgraphs.begin(); p != subgraphs.end(); ++p) {
		if((int)p->size() >= option->min_graph)
			large_graphs.push_back(*p);
		else
			small_graphs.push_back(*p);
	}
	
	// keep only the largest graph if there are only small ones
	if(large_graphs.empty()) {
		auto	max_g = max_element(small_graphs.begin(), small_graphs.end());
		large_graphs.push_back(*max_g);
		small_graphs.erase(max_g);
	}
	
	vector<VCFHeteroHomo *>	subvcfs;
	for(auto p = large_graphs.begin(); p != large_graphs.end(); ++p) {
		subvcfs.push_back(this->make_subvcf(*p));
	}
	
	vector<VCFHeteroHomoRecord *>	unused_records;
	for(auto p = small_graphs.begin(); p != small_graphs.end(); ++p) {
		for(auto q = p->begin(); q != p->end(); ++q)
			unused_records.push_back(this->records[q->first]);
	}
	
	return make_pair(subvcfs, unused_records);
}

string VCFHeteroHomo::make_seq(size_t i) const {
	string	seq;
	for(auto p = records.begin(); p != records.end(); ++p) {
		const auto	*record = *p;
		const int	gt = record->get_which_comes_from(i);
		if(gt == 0)
			seq.push_back('0');
		else if(gt == 1)
			seq.push_back('1');
		else
			seq.push_back('N');
	}
	return seq;
}

bool VCFHeteroHomo::is_all_same_without_N(const string& seq) {
	char	c = '.';
	for(auto p = seq.begin(); p != seq.end(); ++p) {
		if(*p != 'N') {
			if(c == '.')	// initial
				c = *p;
			else if(*p != c)
				return false;
		}
	}
	return true;
}

string VCFHeteroHomo::create_same_color_string(const string& seq) {
	char	c = '0';	// dummy
	for(auto p = seq.begin(); p != seq.end(); ++p) {
		if(*p != 'N') {
			c = *p;
			break;
		}
	}
	return std::string(seq.size(), c);
}

string VCFHeteroHomo::clean_each_sample_seq(size_t i,
								const vector<double>& cMs, double min_c) {
	const string	seq = this->make_seq(i);
	if(is_all_same_without_N(seq))
		return create_same_color_string(seq);
	
	const string	hidden_seq = Imputer::impute(seq, cMs);
	const string	painted_seq = Imputer::paint(hidden_seq, cMs, min_c);
	return painted_seq;
}

void VCFHeteroHomo::clean_each_sample(size_t i, const vector<double>& cMs,
																double min_c) {
	const string	cleaned_seq = clean_each_sample_seq(i, cMs, min_c);
	for(size_t k = 0; k < this->size(); ++k) {
		VCFHeteroHomoRecord	*record = this->records[k];
		const int	w = cleaned_seq.c_str()[k] - '0';
		record->set_int_gt_by_which_comes_from(w, i);
	}
}

void VCFHeteroHomo::clean_in_thread(void *config) {
	auto	*c = (ConfigThreadCleanSeq *)config;
	const size_t	n = c->vcf->num_progenies();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		c->vcf->clean_each_sample(i, c->cMs, c->min_c);
	}
}

void VCFHeteroHomo::clean_each_vcf(const OptionImpute *option) {
	vector<double>	cMs;
	const size_t	L = this->size();
	for(size_t i = 0; i < L; ++i) {
		cMs.push_back(this->record_cM(i));
	}
	
	const double	min_c = option->min_crossover;
	const int	T = option->num_threads;
	vector<ConfigThreadCleanSeq *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThreadCleanSeq(i, T, cMs, min_c, this);
	
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
}

const OptionImpute *VCFHeteroHomo::create_option(int num_threads) const {
	const size_t	num = std::max(2UL, this->size());
	const double	cM_length = this->record_cM(records.size()-1);
	const int	max_dist = std::max(4,
						(int)(cM_length * num_progenies()
										/ num * log10(num) * 2.5 * 0.01));
	return new OptionImpute(max_dist, 20, 5, 1.0, num_threads);
}

pair<vector<VCFHeteroHomo *>, vector<VCFHeteroHomoRecord *>>
										VCFHeteroHomo::clean(int num_threads) {
	if(this->records.empty()) {
		// If it do not create a new VCF,
		// it will delete the VCF that cannot be deleted.
		VCFHeteroHomo	*empty_vcf = new VCFHeteroHomo(header, samples,
								vector<VCFHeteroHomoRecord *>(), get_map());
		return make_pair(vector<VCFHeteroHomo *>(1, empty_vcf),
									vector<VCFHeteroHomoRecord *>());
	}
	
	const OptionImpute	*option = this->create_option(num_threads);
	auto	p = this->determine_haplotype(option);
	vector<VCFHeteroHomo *>&	vcfs = p.first;
	vector<VCFHeteroHomoRecord *>&	unused_records = p.second;
	std::for_each(unused_records.begin(), unused_records.end(),
					std::mem_fun(&VCFImpFamilyRecord::enable_modification));
	for(auto q = vcfs.begin(); q != vcfs.end(); ++q)
		(*q)->clean_each_vcf(option);
	delete option;
	return make_pair(vcfs, unused_records);
}

// How well do the common heteroparents match?
pair<int, int> VCFHeteroHomo::match(const VCFHeteroHomo *other) const {
	if(this->records.empty() || other->records.empty())
		return pair<int, int>(0, 0);
	
	const string	mat_GT1 = this->records.front()->get_GT(0);
	const string	mat_GT2 = other->records.front()->get_GT(0);
	const int	hetero_col1 = mat_GT1.c_str()[0] != mat_GT1.c_str()[2] ? 9 : 10;
	const int	hetero_col2 = mat_GT2.c_str()[0] != mat_GT2.c_str()[2] ? 9 : 10;
	int	num_match = 0;
	int	num_unmatch = 0;
	size_t	k = 0;
	size_t	l = 0;
	while(k < this->size() && l < other->size()) {
		const VCFHeteroHomoRecord	*record1 = this->records[k];
		const VCFHeteroHomoRecord	*record2 = other->records[l];
		if(record1->pos() == record2->pos()) {
			if(record1->get_v()[hetero_col1] == record2->get_v()[hetero_col2])
				num_match += 1;
			else
				num_unmatch += 1;
		}
		if(record1->pos() <= record2->pos())
			k += 1;
		if(record1->pos() >= record2->pos())
			l += 1;
	}
	
	return make_pair(num_match, num_unmatch);
}

void VCFHeteroHomo::inverse_hetero_parent_phases() {
	const int	hetero_index = this->is_mat_hetero() ? 0 : 1;
	for(auto p = this->records.begin(); p != this->records.end(); ++p) {
		if((*p)->get_GT(hetero_index) == "0|1")
			(*p)->set_GT(hetero_index, "1|0");
		else
			(*p)->set_GT(hetero_index, "0|1");
	}
}

// Create a pair of VCFHeteroHomo for each Family and store it for each parent
tuple<VCFHeteroHomo *, VCFHeteroHomo *, vector<VCFHeteroHomoRecord *>>
VCFHeteroHomo::make_VCFHeteroHomo(const vector<VCFHeteroHomoRecord *>& records,
								  const vector<STRVEC>& header,
								  const STRVEC& samples, const Map& geno_map) {
	// classify records
	vector<VCFHeteroHomoRecord *>	heho_mat_records;
	vector<VCFHeteroHomoRecord *>	heho_pat_records;
	vector<VCFHeteroHomoRecord *>	unused_records;
	for(auto p = records.begin(); p != records.end(); ++p) {
		VCFHeteroHomoRecord	*record = *p;
		if(record == NULL)
			continue;
		else if(!record->is_imputable())
			unused_records.push_back(record);
		else if(record->is_mat_hetero())
			heho_mat_records.push_back(record);
		else
			heho_pat_records.push_back(record);
	}
	
	auto	*vcf_mat = new VCFHeteroHomo(header, samples,
											heho_mat_records, geno_map);
	auto	*vcf_pat = new VCFHeteroHomo(header, samples,
											heho_pat_records, geno_map);
	return make_tuple(vcf_mat, vcf_pat, unused_records);
}

// collect VCFs whose hetero parent is the same, and impute it
// and change the phase to be as same as possible
ImpResult VCFHeteroHomo::clean_vcfs(
							const vector<VCFHeteroHomoRecord *>& records,
							const vector<STRVEC>& header, const STRVEC& samples,
							const Map& geno_map, int num_threads) {
	auto	tuple1 = make_VCFHeteroHomo(records, header, samples, geno_map);
	auto	*vcf_mat = get<0>(tuple1);
	auto	*vcf_pat = get<1>(tuple1);
	auto	unused = get<2>(tuple1);
	
	auto	pair1 = vcf_mat->clean(num_threads);
	auto	pair2 = vcf_pat->clean(num_threads);
	auto	vcfs_mat = pair1.first;
	auto	unused_records1 = pair1.second;
	auto	vcfs_pat = pair2.first;
	auto	unused_records2 = pair2.second;
	vcf_mat->clear_records();
	delete vcf_mat;
	vcf_pat->clear_records();
	delete vcf_pat;
	
	vector<VCFHeteroHomo *>	vcfs;
	vcfs.insert(vcfs.end(), vcfs_mat.begin(), vcfs_mat.end());
	vcfs.insert(vcfs.end(), vcfs_pat.begin(), vcfs_pat.end());
	unused.insert(unused.end(), unused_records1.begin(), unused_records1.end());
	unused.insert(unused.end(), unused_records2.begin(), unused_records2.end());
	return make_pair(vcfs, unused);
}

void VCFHeteroHomo::inverse_phases(const vector<VCFHeteroHomo *>& vcfs) {
	
	const auto	*graph = make_vcf_graph(vcfs);
	const vector<bool>	invs = graph->optimize_inversions();
	delete graph;
	for(size_t i = 0; i < vcfs.size(); ++i) {
		if(invs[i])
			vcfs[i]->inverse_hetero_parent_phases();
	}
}

const InverseGraph *VCFHeteroHomo::make_vcf_graph(
									const vector<VCFHeteroHomo *>& vcfs) {
	InverseGraph	*graph = new InverseGraph();
	const size_t	N = vcfs.size();
	for(size_t i = 0; i < N; ++i) {
		(*graph)[i];	// initialize vector
		for(size_t j = i + 1; j < N; ++j) {
			const auto	p = vcfs[i]->match(vcfs[j]);
			const int	num_match = p.first;
			const int	num_unmatch = p.second;
			if(num_match != 0 || num_unmatch != 0) {
				(*graph)[i].push_back(make_tuple(j, num_match, num_unmatch));
				(*graph)[j].push_back(make_tuple(i, num_match, num_unmatch));
			}
		}
	}
	return graph;
}
