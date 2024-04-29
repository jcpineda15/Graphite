#ifndef __IMPUTER
#define __IMPUTER

#include <memory>
#include "Baum_Welch_with_fixed_Ts.h"


namespace Imputer {

using Color = int;


//////////////////// List ////////////////////

// able to add an item and immutable
template<typename T>
class List {
public:
	using ptrList = std::shared_ptr<const List<T>>;
	
private:
	ptrList	first;
	T		last;
	
public:
	List<T>(T val) : first(NULL), last(val) { }
	List<T>(const ptrList f, T& val) : first(f), last(val) { }
	~List<T>() { }
	
	const ptrList prev() const { return first; }
	T back() const { return last; }
	bool is_null() const { return first == NULL; }
	std::vector<T> to_vector() const;
	
public:
	static ptrList add(const ptrList ptr, T val);
};


//////////////////// Region ////////////////////

class Region {
	const Color			color;
	const std::size_t	first;
	const std::size_t	last;
	const double		cM;
	
public:
	Region(Color c, std::size_t f, std::size_t l, double cm) :
							color(c), first(f), last(l), cM(cm) { }
	
	Color get_color() const { return color; }
	std::size_t get_first() const { return first; }
	std::size_t get_last() const { return last; }
	double get_cM() const { return cM; }
	std::size_t get_size() const { return last - first; }
	std::string str() const;
	
public:
	static std::vector<const Region *> create(const std::string& seq,
											const std::vector<double>& cMs);
	static std::vector<std::tuple<Color, std::size_t, std::size_t>>
										group_by_colors(const std::string& seq);
};


//////////////////// State ////////////////////

class State;
using ptrState = std::shared_ptr<const State>;
using ptrColors = List<Color>::ptrList;

class State {
	const std::vector<const Region *>&	regions;
	ptrColors	color_list;
	const std::size_t	size;
	const int			num_continuous;
	const double		continuous_length;
	const double		painted_len;
	
public:
	State(const std::vector<const Region *>& rs,
			ptrColors cs, std::size_t s,
			int nc, double cl, double pl) :
				regions(rs), color_list(cs), size(s), num_continuous(nc),
				continuous_length(cl), painted_len(pl) { }
	State(const std::vector<const Region *>& rs);
								
	~State() { }
	
	int get_num_continuous() const { return num_continuous; }
	double painted_length() const { return painted_len; }
	Color last_color() const { return color_list->back(); }
	Color current_color() const { return regions[size-1]->get_color(); }
	Color next_color() const { return regions[size]->get_color(); }
	bool is_extendable() const {
		return size != regions.size() - 1 || next_color() ==  last_color();
	}
	bool is_changable(double MIN_CROSSOVER) const;
	std::vector<Color> get_colors() const;
	std::string str() const;
	
	State *add(double cM) const;
	State *change(double cM) const;
	std::vector<const State *> nexts(double cM, double MIN_CROSSOVER) const;
	
public:
	static ptrState min(const std::vector<ptrState>& states);
	static std::vector<ptrState> select(
						const std::vector<const State *>& states);
	static std::vector<ptrState> next(std::vector<ptrState>& states,
										double cM, double MIN_CROSSOVER);
};


//////////////////// Paint ////////////////////

std::string paint(const std::string& seq, const std::vector<double>& cMs,
													double MIN_CROSSOVER);


//////////////////// impute ////////////////////

std::vector<char> create_states(const std::string& seq);
BaumWelch::TransitionMatrix compute_T(double p);
std::string impute(const std::string& seq, const std::vector<double>& cMs);

}
#endif
