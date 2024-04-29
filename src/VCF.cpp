#include <algorithm>
#include "../include/VCF.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFRecord ////////////////////

pair<string,ll> VCFRecord::position() const {
	return pair<string,ll>(v[0], stoll(v[1]));
}

string VCFRecord::gt(const string& sample) const {
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		if(*p == sample)
			return v[p-samples.begin()+9];
	}
	return "";
}

STRVEC VCFRecord::gts() const {
	STRVEC	w(samples.size());
	std::copy(v.begin() + 9, v.end(), w.begin());
	return w;
}

#include <cassert>
vector<int> VCFRecord::get_int_gts() const {
	vector<int>	w(samples.size());
	for(size_t i = 0; i < samples.size(); ++i) {
		assert(i < w.size());
		w[i] = get_int_gt(i);
	}
	return w;
}

bool VCFRecord::is_homo(size_t i) const {
	const string&	s = this->v[i+9];
	const char	*ptr = s.c_str();
	return ptr[0] != '.' && ptr[0] == ptr[2];
}

bool VCFRecord::is_hetero(size_t i) const {
	const string&	s = this->v[i+9];
	const char	*ptr = s.c_str();
	return ptr[0] != '.' && ptr[2] != '.' && ptr[0] != ptr[2];
}

STRVEC VCFRecord::extract_v(const STRVEC& samples) const {
	map<string,size_t>	dic;
	for(size_t i = 0U; i < this->samples.size(); ++i)
		dic[this->samples[i]] = i;
	
	STRVEC	new_v(v.begin(), v.begin() + 9);
	for(auto p = samples.begin(); p != samples.end(); ++p)
		new_v.push_back(v[dic[*p]+9]);
	return new_v;
}

void VCFRecord::write(ostream& os) const {
	Common::write_tsv(this->v, os);
}

void VCFRecord::copy_properties(STRVEC::iterator it) const {
	std::copy(v.begin(), v.begin() + 9, it);
}

void VCFRecord::set_GT(size_t i, const string& gt) {
	this->v[i+9].replace(0, 3, gt);
}

void VCFRecord::set_int_GT(size_t i, int gt) {
	switch(gt) {
		case  0: set_GT(i, "0/0"); return;
		case  1: set_GT(i, "0/1"); return;
		case  2: set_GT(i, "1/1"); return;
		default: set_GT(i, "./."); return;
	}
}


//////////////////// VCFBase ////////////////////

VCFBase::VCFBase(const vector<STRVEC>& h, const STRVEC& s) : 
								header(h), samples(s),
								sample_ids(number_samples(s)) {
	determine_chromosome_id();
}

map<string,size_t> VCFBase::number_samples(const STRVEC& samples_) const {
	map<string,size_t>	dic;
	size_t	id = 0U;
	for(auto p = samples_.begin(); p != samples_.end(); ++p) {
		dic[*p] = id;
		++id;
	}
	return dic;
}

POSITION VCFBase::position(const pair<string,ll>& p) const {
	const string&	chr = p.first;
	const ll		pos = p.second;
	auto	q = this->chrs.find(chr);
	if(q == this->chrs.end()) {
		this->chrs[chr] = (int)this->chrs.size() + 1;
	}
	return make_pair(this->chrs[chr], pos);
}

POSITION VCFBase::record_position(const VCFRecord& record) const {
	return this->position(record.position());
}

std::string VCFBase::chr(int chr_id) const {
	for(auto p = this->chrs.begin(); p != this->chrs.end(); ++p) {
		if(p->second == chr_id)
			return p->first;
	}
	return "";
}

int VCFBase::find_column(const string& sample) const {
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		if(*p == sample)
			return (p - samples.begin()) + 9;
	}
	return -1;
}

void VCFBase::write_header(ostream& os) const {
	for(auto p = header.begin(); p != header.end(); ++p)
		Common::write_tsv(*p, os);
}

void VCFBase::determine_chromosome_id() {
	int	id = 0;
	for(auto p = header.begin(); p != header.end(); ++p) {
		const string&	s = p->front();
		if(s.find("##config=<ID=") == string::npos)
			continue;
		size_t	q = s.find(",length");
		string	chr	= s.substr(13U, q);
		chrs[chr] = id;
		id += 1;
	}
}


//////////////////// VCFReader ////////////////////

VCFReader::VCFReader(const string& path) {
	if(path.substr(path.length() - 3U) != ".gz")
		reader = new FileReader(path);
	else
		reader = new FileReaderGZ(path);
}

VCFReader::~VCFReader() {
	delete reader;
}

void VCFReader::read_header() {
	string	line;
	while(reader->getline(line)) {
		header.push_back(Common::split(line, '\t'));
		if(line.substr(0U, 6U) == "#CHROM")
			break;
	}
}

STRVEC VCFReader::next() {
	string	line;
	if(!reader->getline(line))
		return STRVEC();
	
	return Common::split(line, '\t');
}

STRVEC VCFReader::get_samples() const {
	const STRVEC&	v = header.back();
	STRVEC	samples(v.begin() + 9, v.end());
	return samples;
}


//////////////////// VCFSmallBase ////////////////////

vector<STRVEC> VCFSmallBase::trim_header(const STRVEC& samples) const {
	vector<STRVEC>	header = this->get_header();
	STRVEC	v(header.back().begin(), header.back().begin() + 9);
	v.insert(v.end(), samples.begin(), samples.end());
	header.back() = v;
	return header;
}

vector<int> VCFSmallBase::clip_raw_haplotype(size_t sample_id, int i) const {
	vector<int>	hap;
	for(size_t j = 0; j < size(); ++j) {
		const VCFRecord	*record = get_record(j);
		const string&	gt = record->get_gt(sample_id);
		const char	c = gt.c_str()[i*2];
		hap.push_back(c == '.' ? -1 : (int)(c - '0'));
	}
	return hap;
}

vector<size_t> VCFSmallBase::extract_columns(STRVEC::const_iterator first,
											STRVEC::const_iterator last) const {
	const STRVEC&	samples = this->get_samples();
	map<string, size_t>	dic;
	for(size_t i = 0; i < samples.size(); ++i)
		dic[samples[i]] = i + 9;
	
	vector<size_t>	columns;
	for(auto p = first; p != last; ++p) {
		auto	q = dic.find(*p);
		if(q != dic.end())
			columns.push_back(q->second);
		else
			columns.push_back(string::npos);
	}
	return columns;
}

VCFSmall *VCFSmallBase::extract_samples(const STRVEC& samples) const {
	const auto	header = trim_header(samples);
	vector<VCFRecord *>	empty_records;
	VCFSmall	*vcf = new VCFSmall(header, samples, empty_records);
	const vector<size_t>	cs = extract_columns(samples);
	for(size_t i = 0; i < size(); ++i) {
		const VCFRecord	*record = get_record(i);
		const STRVEC&	v = record->get_v();
		STRVEC	new_v(v.begin(), v.begin() + 9);
		for(auto q = cs.begin(); q != cs.end(); ++q)
			new_v.push_back(v[*q]);
		VCFRecord	*new_record = new VCFRecord(new_v, vcf->get_samples());
		vcf->add_record(new_record);
	}
	return vcf;
}

void VCFSmallBase::write(ostream& os, bool write_header) const {
	if(write_header)
		this->write_header(os);
	for(size_t i = 0; i < size(); ++i)
		get_record(i)->write(os);
}

void VCFSmallBase::write_header(ostream& os) const {
	const auto&	header = get_header();
	for(auto p = header.begin(); p != header.end(); ++p)
		Common::write_tsv(*p, os);
}


//////////////////// VCFSmall ////////////////////

VCFSmall::VCFSmall(const vector<STRVEC>& h, const STRVEC& s,
										vector<VCFRecord *> rs, bool rr) :
				VCFBase(h, s), VCFSmallBase(), records(rs), reuses_records(rr) {
	for(auto p = records.begin(); p != records.end(); ++p)
		this->record_position(**p);
}

VCFSmall::~VCFSmall() {
	if(!reuses_records) {
		for(auto p = records.begin(); p != records.end(); ++p)
			delete *p;
	}
}

VCFSmall *VCFSmall::read(const string& path) {
	VCFReader	*reader = new VCFReader(path);
	reader->read_header();
	const vector<STRVEC>	header = reader->get_header();
	const STRVEC	samples = reader->get_samples();
	
	vector<VCFRecord *>	records;
	while(true) {
		const STRVEC	v = reader->next();
		if(v.empty())
			break;
		VCFRecord	*record = new VCFRecord(v, samples);
		records.push_back(record);
	}
	delete reader;
	
	return new VCFSmall(header, samples, records);
}

// join VCFs in order of given samples
VCFSmall *VCFSmall::join(const vector<const VCFSmallBase *>& vcfs,
													const STRVEC& samples) {
	map<string, pair<const VCFSmallBase *, size_t>>	dic;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const VCFSmallBase	*vcf = *p;
		const STRVEC&	ss = vcf->get_samples();
		for(size_t i = 0; i < ss.size(); ++i)
			dic[ss[i]] = make_pair(vcf, i + 9);
	}
	
	vector<tuple<string, const VCFSmallBase *, size_t>>	cols;
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		auto	q = dic.find(*p);
		if(q != dic.end())
			cols.push_back(make_tuple(*p, q->second.first, q->second.second));
	}
	
	STRVEC	new_samples;
	for(auto p = cols.begin(); p != cols.end(); ++p)
		new_samples.push_back(get<0>(*p));
	
	// new_vcf owns samples
	auto	new_header = vcfs.front()->trim_header(new_samples);
	vector<VCFRecord *>	empty_records;
	VCFSmall	*new_vcf = new VCFSmall(new_header, new_samples, empty_records);
	const STRVEC&	samples_ = new_vcf->get_samples();
	for(size_t i = 0; i < vcfs.front()->size(); ++i) {
		const VCFSmallBase	*vcf = vcfs[0];
		const STRVEC&	orig_v = vcf->get_record(i)->get_v();
		STRVEC	v(orig_v.begin(), orig_v.begin() + 9);
		for(auto p = cols.begin(); p != cols.end(); ++p)
			v.push_back(get<1>(*p)->get_record(i)->get_v()[get<2>(*p)]);
		new_vcf->add_record(new VCFRecord(v, samples_));
	}
	
	return new_vcf;
}

VCFSmall *VCFSmall::join(const VCFSmallBase *vcf1, const VCFSmallBase *vcf2,
														const STRVEC& samples) {
	vector<const VCFSmallBase *>	vcfs{ vcf1, vcf2 };
	return VCFSmall::join(vcfs, samples);
}


//////////////////// VCFHuge::ChromDivisor ////////////////////

VCFSmall *VCFHuge::ChromDivisor::next() {
	while(true) {
		VCFRecord	*record = this->vcf->next();
		if(state == STATE::START) {
			if(record == NULL)
				return NULL;
			else {
				this->chrom = record->chrom();
				this->records.push_back(record);
				state = STATE::DOING;
			}
		}
		else if(state == STATE::END) {
			return NULL;
		}
		else if(record == NULL) {
			VCFSmall	*vcf_chrom = new VCFSmall(this->vcf->header,
											this->vcf->samples, this->records);
			this->records.clear();
			state = STATE::END;
			return vcf_chrom;
		}
		else if(record->chrom() != this->chrom) {
			VCFSmall	*vcf_chrom = new VCFSmall(this->vcf->header,
											this->vcf->samples, this->records);
			this->records.clear();
			this->records.push_back(record);
			this->chrom = record->chrom();
			return vcf_chrom;
		}
		else {
			this->records.push_back(record);
		}
	}
}


//////////////////// VCFHuge ////////////////////

VCFHuge::VCFHuge(const vector<STRVEC>& h, const STRVEC& s, VCFReader *r_) :
												VCFBase(h, s), reader(r_) { }

VCFHuge::~VCFHuge() {
	delete reader;
}

VCFRecord *VCFHuge::next() {
	const STRVEC	v = reader->next();
	if(v.empty())
		return NULL;
	
	VCFRecord	*record = new VCFRecord(v, this->samples);
	this->record_position(*record);		// necessary to update chrs
	return record;
}

VCFRecord *VCFHuge::proceed(POSITION pos) {
	while(true) {
		VCFRecord	*record = next();
		if(record == NULL)
			return NULL;
		
		if(this->record_position(*record) == pos)
			return record;
		else
			delete record;
	}
	return NULL;
}

VCFHuge *VCFHuge::read(const string& path) {
	VCFReader	*reader = new VCFReader(path);
	reader->read_header();
	const vector<STRVEC>&	header = reader->get_header();
	const STRVEC	samples = reader->get_samples();
	return new VCFHuge(header, samples, reader);
}
