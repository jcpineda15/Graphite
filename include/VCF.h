#ifndef __VCF
#define __VCF

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "filereader.h"
#include "Genotype.h"


typedef long long					ll;
typedef std::vector<std::string>	STRVEC;
typedef std::pair<int, ll>			POSITION;


//////////////////// VCFRecord ////////////////////

class VCFRecord {
protected:
	STRVEC	v;
	const STRVEC&	samples;
	
public:
	VCFRecord(const STRVEC& v_, const STRVEC& s_) : v(v_), samples(s_) { }
	virtual ~VCFRecord() { }
	
	const STRVEC& get_v() const { return v; }
	const STRVEC& get_samples() const { return samples; }
	size_t num_samples() const { return samples.size(); }
	const std::string& chrom() const { return this->v[0]; }
	std::pair<std::string,ll>	position() const;
	ll	pos() const { return stoll(this->v[1]); }
	const std::string	format() const { return this->v[8]; }
	std::string	gt(const std::string& sample) const;
	bool is_NA(std::size_t i) const { return Genotype::is_NA(v[i+9]); }
	STRVEC gts() const;
	const std::string& get_gt(std::size_t i) const { return v[i+9]; }
	std::string& get_mut_gt(std::size_t i) { return v[i+9]; }
	std::string get_GT(std::size_t i) const { return v[i+9].substr(0, 3); }
	int get_int_gt(std::size_t i) const { return Genotype::get_int_gt(v[i+9]); }
	std::vector<int> get_int_gts() const;
	bool is_homo(std::size_t i) const;
	bool is_hetero(std::size_t i) const;
	STRVEC extract_v(const STRVEC& samples) const;
	void write(std::ostream& os) const;
	
	void copy_properties(STRVEC::iterator it) const;
	void set_GT(std::size_t i, const std::string& gt);
	void set_int_GT(std::size_t i, int gt);
	void set(const STRVEC& new_v) { v = new_v; }
};


//////////////////// VCFBase ////////////////////

class VCFBase {
protected:
	const std::vector<STRVEC>	header;
	const STRVEC	samples;
	const std::map<std::string,std::size_t>	sample_ids;
	mutable std::map<std::string,int>	chrs;
	
public:
	VCFBase(const std::vector<STRVEC>& header_, const STRVEC& samples_);
	virtual ~VCFBase() { }
	
	const std::vector<STRVEC>& get_header() const { return header; }
	const STRVEC& get_samples() const { return samples; }
	size_t num_samples() const { return samples.size(); }
	std::map<std::string,std::size_t>
	number_samples(const STRVEC& samples_) const;
	POSITION position(const std::pair<std::string,ll>& p) const;
	POSITION record_position(const VCFRecord& record) const;
	std::string chr(int chr_id) const;
	int find_column(const std::string& sample) const;
	void write_header(std::ostream& os) const;
	
	void copy_chrs(VCFBase *vcf) const { vcf->chrs = chrs; }
	
private:
	void determine_chromosome_id();
};


//////////////////// VCFReader ////////////////////

class VCFReader {
	FileReaderBase	*reader;
	std::vector<STRVEC>	header;
	
public:
	VCFReader(const std::string& path);
	~VCFReader();
	
	void read_header();
	const std::vector<STRVEC>& get_header() { return header; }
	STRVEC get_samples() const;
	STRVEC next();
};


//////////////////// VCFSmallBase ////////////////////

class VCFSmall;

class VCFSmallBase {
public:
	VCFSmallBase() { }
	virtual ~VCFSmallBase() { }
	
	virtual const std::vector<STRVEC>& get_header() const = 0;
	virtual const STRVEC& get_samples() const = 0;
	virtual std::size_t size() const = 0;
	virtual VCFRecord *get_record(std::size_t i) const = 0;
	
	std::vector<STRVEC> trim_header(const STRVEC& samples) const;
	std::vector<int> clip_raw_haplotype(std::size_t sample_id, int i) const;
	std::vector<std::size_t> extract_columns(const STRVEC& samples) const {
		return extract_columns(samples.begin(), samples.end());
	}
	std::vector<std::size_t> extract_columns(STRVEC::const_iterator first,
											STRVEC::const_iterator last) const;
	VCFSmall *extract_samples(const STRVEC& samples) const;
	void write(std::ostream& os, bool write_header=true) const;
	void write_header(std::ostream& os) const;
};


//////////////////// VCFSmall ////////////////////

class VCFSmall : public VCFBase, public VCFSmallBase {
protected:
	std::vector<VCFRecord *>	records;
	bool	reuses_records;
	
public:
	VCFSmall(const std::vector<STRVEC>& header, const STRVEC& samples,
								std::vector<VCFRecord *> rs, bool rr=false);
	virtual ~VCFSmall();
	
	///// virtual methods /////
	const std::vector<STRVEC>& get_header() const {
		return VCFBase::get_header();
	}
	const STRVEC& get_samples() const { return VCFBase::get_samples(); }
	std::size_t size() const { return records.size(); }
	VCFRecord *get_record(std::size_t i) const { return records[i]; }
	
	///// non-virtual methods /////
	const std::vector<VCFRecord *>& get_records() const { return records; }
	bool is_empty() const { return records.empty(); }
	
	void add_record(VCFRecord *record) { records.push_back(record); }
	void add_records(std::vector<VCFRecord *>& rs) {
		records.insert(records.end(), rs.begin(), rs.end());
	}
	void clear_records() { records.clear(); }
	
public:
	static VCFSmall *read(const std::string& path);
	// join VCFs in order of given samples
	static VCFSmall *join(const std::vector<const VCFSmallBase *>& vcfs,
														const STRVEC& samples);
	static VCFSmall *join(const VCFSmallBase *vcf1, const VCFSmallBase *vcf2,
														const STRVEC& samples);
};


//////////////////// VCFHuge ////////////////////

class VCFHuge : public VCFBase {
public:
	class ChromDivisor {
		enum class STATE { START, DOING, END };
		VCFHuge	*vcf;
		std::vector<VCFRecord *>	records;
		std::string	chrom;
		STATE	state;
		
	public:
		ChromDivisor(VCFHuge *v) : vcf(v), state(STATE::START) { }
		VCFSmall *next();
	};
	
private:
	VCFReader	*reader;
	
public:
	VCFHuge(const std::vector<STRVEC>& header_, const STRVEC& samples_,
													VCFReader *reader_);
	~VCFHuge();
	
	VCFRecord *next();
	VCFRecord *proceed(POSITION pos);
	
public:
	static VCFHuge *read(const std::string& path);
};
#endif
