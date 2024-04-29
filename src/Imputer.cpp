#include <map>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "../include/Imputer.h"
#include "../include/common.h"

using namespace std;
using namespace Imputer;


//////////////////// List ////////////////////

template<typename T>
vector<T> List<T>::to_vector() const {
	vector<T>	v;
	const List<T>	*ptr_list = this;
	while(ptr_list != NULL) {
		v.push_back(ptr_list->last);
		ptr_list = ptr_list->first.get();
	}
	
	std::reverse(v.begin(), v.end());
	return v;
}

template<typename T>
typename List<T>::ptrList List<T>::add(
											const ptrList ptr, T val) {
	return ptrList(new List<T>(ptr, val));
}


//////////////////// Region ////////////////////

string Region::str() const {
	stringstream	ss;
	ss << '[' << first << ", " << last << ']';
	return ss.str();
}

vector<tuple<Color, size_t, size_t>> Region::group_by_colors(
														const string& seq) {
	size_t	last = 0U;
	Color	prev_color = -1;
	vector<tuple<Color, size_t, size_t>>	intervals;	// [(Color, first, last)]
	for(auto p = seq.begin(); p != seq.end(); ++p) {
		const Color	color = (Color)(*p) - '0';
		if(prev_color != -1 && color != prev_color) {
			const size_t	j = p - seq.begin();
			intervals.push_back(make_tuple(prev_color, last, j));
			last = p - seq.begin();
		}
		prev_color = color;
	}
	intervals.push_back(make_tuple(prev_color, last, seq.size()));
	return intervals;
}

vector<const Region *> Region::create(const string& seq,
											const vector<double>& cMs) {
	const auto	intervals = group_by_colors(seq);
	
	vector<const Region *>	regions;
	for(auto p = intervals.begin(); p != intervals.end(); ++p) {
		const Color		color = get<0>(*p);
		const size_t	first = get<1>(*p);
		const size_t	last = get<2>(*p);
		const double	cM1 = first == 0 ? cMs[first] :
											(cMs[first-1] + cMs[first]) / 2;
		const double	cM2 = last == seq.size() ? cMs[last-1] :
											(cMs[last-1] + cMs[last]) / 2;
		Region	*region = new Region(color, first, last, cM2 - cM1);
		regions.push_back(region);
	}
	return regions;
}


//////////////////// State ////////////////////

State::State(const std::vector<const Region *>& rs) :
						regions(rs),
						color_list(new List<Color>(regions[0]->get_color())),
						size(1U), num_continuous(1U),
						continuous_length(regions[0]->get_cM()),
						painted_len(0.0) { }

bool State::is_changable(double MIN_CROSSOVER) const {
	if(continuous_length >= MIN_CROSSOVER)
		return true;
	
	// it is painted the same color till now
	ptrColors	colors = color_list;
	while(true) {
		if(colors->back() != last_color())
			return false;
		if(colors->is_null())
			break;
		colors = colors->prev();
	}
	return true;
}

State *State::add(double cM) const {
	const bool	painted = last_color() != next_color();
	ptrColors	colors = List<Color>::add(color_list, last_color());
	const double	new_painted_length = painted_len + (painted ? cM : 0.0);
	return new State(regions, colors, size + 1, num_continuous + 1,
							continuous_length + cM, new_painted_length);
}

State *State::change(double cM) const {
	const Color	color = last_color() == 0 ? 1 : 0;
	const bool	painted = color != next_color();
	ptrColors	colors = List<Color>::add(color_list, color);
	const double	new_painted_length = painted_len + (painted ? cM : 0.0);
	return new State(regions, colors, size + 1, 1, cM, new_painted_length);
}

vector<Color> State::get_colors() const {
	vector<Color>	colors;
	vector<Color>	region_colors = color_list.get()->to_vector();
	for(auto p = region_colors.begin(); p != region_colors.end(); ++p) {
		const size_t	k = p - region_colors.begin();
		const Region	*region = *(regions.begin() + k);
		for(size_t i = 0U; i < region->get_size(); ++i)
			colors.push_back(*p);
	}
	return colors;
}

vector<const State *> State::nexts(double cM, double MIN_CROSSOVER) const {
	vector<const State *>	states;
	if(is_extendable())
		states.push_back(add(cM));
	if(is_changable(MIN_CROSSOVER))
		states.push_back(change(cM));
	return states;
}

ptrState State::min(const vector<ptrState>& states) {
	ptrState	min_state = states.front();
	for(auto p = states.begin() + 1; p != states.end(); ++p) {
		const State	*state = p->get();
		if(state->painted_length() < min_state.get()->painted_length())
			min_state = *p;
	}
	return min_state;
}

// minimize states
vector<ptrState> State::select(const vector<const State *>& states) {
	typedef pair<int,Color>	Key;
	map<Key,vector<ptrState>>	dic;
	for(auto p = states.begin(); p != states.end(); ++p) {
		Key	key((*p)->get_num_continuous(),(*p)->last_color());
		dic[key].push_back(ptrState(*p));
	}
	
	vector<ptrState>	selected_states;
	for(auto p = dic.begin(); p != dic.end(); ++p)
		selected_states.push_back(min(p->second));
	return selected_states;
}

vector<ptrState> State::next(vector<ptrState>& states,
										double cM, double MIN_CROSSOVER) {
	vector<const State *>	new_states;
	for(auto p = states.begin(); p != states.end(); ++p) {
		const auto	next_states = p->get()->nexts(cM, MIN_CROSSOVER);
		new_states.insert(new_states.end(),
							next_states.begin(), next_states.end());
	}
	return select(new_states);
}

string State::str() const {
	stringstream	ss;
	for(auto p = regions.begin(); p != regions.end(); ++p)
		ss << (*p)->str();
	ss << "\n";
	
	const vector<Color>	colors = this->get_colors();
	for(auto p = colors.begin(); p != colors.end(); ++p)
		ss << *p;
	
	return ss.str();
}


//////////////////// Paint ////////////////////

string Imputer::paint(const string& seq, const vector<double>& cMs,
												double MIN_CROSSOVER) {
	const auto	regions = Region::create(seq, cMs);
	if(regions.size() <= 2U) {
		Common::delete_all(regions);
		return seq;
	}
	
	auto	states = vector<ptrState>(1U, ptrState(new State(regions)));
	for(auto p = regions.begin() + 1; p != regions.end(); ++p)
		states = State::next(states, (*p)->get_cM(), MIN_CROSSOVER);
	
	ptrState	state = State::min(states);
	const auto	colors = state.get()->get_colors();
	
	stringstream	ss;
	for(auto p = colors.begin(); p != colors.end(); ++p)
		ss << *p;
	Common::delete_all(regions);
	return ss.str();
}


//////////////////// impute ////////////////////

vector<char> Imputer::create_states(const string& seq) {
	if(seq.find('N') != string::npos) {
		vector<char>	states = { '0', '1', 'N' };
		return states;
	}
	else {
		vector<char>	states = { '0', '1' };
		return states;
	}
}

BaumWelch::TransitionMatrix Imputer::compute_T(double prob) {
	BaumWelch::TransitionMatrix	T;
	for(size_t h0 = 0; h0 < 2; ++h0) {
		for(size_t h1 = 0; h1 < 2; ++h1)
			T[h0*2+h1] = h0 == h1 ? 1.0 - prob : prob;
	}
	return T;
}

string Imputer::impute(const string& seq, const vector<double>& cMs) {
	// Kosambi
	vector<double>	ps;
	for(size_t i = 0; i < cMs.size() - 1; ++i) {
		const double	d = (cMs[i+1] - cMs[i]) / 100;
		const double	r = (exp(d) - exp(-d)) / (exp(d) + exp(-d));
		ps.push_back(r);
	}
	
	vector<BaumWelch::TransitionMatrix>	Ts;
	for(auto p = ps.begin(); p != ps.end(); ++p)
		Ts.push_back(Imputer::compute_T(*p));
	
	const string	hidden_seq = BaumWelch::impute(seq, Ts);
	return hidden_seq;
}
