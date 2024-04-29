# coding: utf-8
# Baum-Welch_with_fixed_Ts_log.py
# 1次元の列をimputeする

from functools import reduce
from itertools import count, product
from collections import defaultdict
from math import log, exp, floor
from Matrix import Matrix


#################### Log ####################

LOGZERO = -10000.0

def log_add(log_a, log_b):
	# e.g. log_a = log2 log_b = log3
	# log(2 + 3) = log3 + log(2/3 + 1) = log3 + log(exp(log2-log3) + 1)
	if log_a < log_b:
		return log_b + log(exp(log_a - log_b) + 1.0)
	else:
		return log_a + log(exp(log_b - log_a) + 1.0)

def log_sub(log_a, log_b):
	# e.g. log_a = log3 log_b = log2
	# log(3 - 2) = log3 + log(1 - 2/3) = log3 + log(1 - exp(log2-log3))
	assert(log_a >= log_b)
	if log_a == log_b or exp(log_b - log_a) == 1.0:
		return LOGZERO
	
	return log_a + log(1.0 - exp(log_b - log_a))

def log_sum(iterable):
	return reduce(log_add, iterable)

def log_diff(log_a, log_b):
	if log_a > log_b:
		return log_sub(log_a, log_b)
	else:
		return log_sub(log_b, log_a)

def mat_log(M):
	return dict((k, log(v) if v > 0.0 else LOGZERO) for k, v in M.items())

def modified_log(x):
	return -10000.0 if x == 0.0 else log(x)


#################### Matrix ####################

# M**e e: float
def power(M, log_M, e):
	if e < 1:
		return (log_M * e).exp()
	
	e1 = int(floor(e))
	return M**e1 * (log_M * (e - e1)).exp()


#################### Baum-Welch ####################

def initialize_pi(hidden_states):
	prob = -log(len(hidden_states))
	return dict((h, prob) for h in hidden_states)

def initialize_transition_matrix(hidden_states):
	def prob(h1, h2):
		if h1 == h2:
			return log(0.9)
		else:
			return log(0.1 / (N - 1))
	
	N = len(hidden_states)
	return dict(((h1, h2), prob(h1, h2))
					for h1, h2 in product(hidden_states, repeat=2))

def initialize_emission_matrix(hidden_states, states):
	def prob(h, e):
		if h == e:
			return log(0.9)
		else:
			return log(0.1 / (N - 1))
	
	N = len(states)
	return dict(((h, e), prob(h, e)) for h, e in product(hidden_states, states))

def compute_parameters(seq, hidden_states, pi, Ts, A):
	L = len(seq)
	N = len(hidden_states)
	
	# compute α
	a = [ defaultdict(float) for s in seq ]		# probability of s
	for h in hidden_states:
		a[0][h] = pi[h]+A[(seq[0],h)]
	for t, s in enumerate(seq[1:], 1):
		T = Ts[t-1]
		for h in hidden_states:
			a[t][h] = log_sum(a[t-1][h_prev]+T[(h_prev,h)]+A[(s,h)]
										for h_prev in hidden_states)
	
	# compute β
	b = [ defaultdict(float) for s in seq ]
	for h in hidden_states:
		b[-1][h] = 0.0	# log(1.0)
	for t, s in zip(range(L-2, -1, -1), seq[-1::-1]):
		T = Ts[t]
		for h in hidden_states:
			b[t][h] = log_sum(T[(h,h_next)]+A[(s,h_next)]+b[t+1][h_next]
										for h_next in hidden_states)
	
	# compute γ
	gamma1 = [ dict((h, a1[h]+b1[h]) for h in hidden_states)
										for a1, b1 in zip(a, b) ]
	prob_seq = log_sum(gamma1[0].values())	# probability that seq is observed
	gamma = [ dict((h, p-prob_seq) for h, p in g.items()) for g in gamma1 ]
	
	# compute ξ
	xi = [ defaultdict(float) for s in seq[:-1] ]
	for t, s in enumerate(seq[1:]):
		T = Ts[t]
		for hi, hj in product(hidden_states, repeat=2):
			xi[t][(hi,hj)] = a[t][hi]+T[(hi,hj)]+A[(s,hj)]+b[t+1][hj] - prob_seq
	
	return (a, b, xi, gamma)

def reverse_probs(A, hidden_states, states):
	sum_observed = dict((s0, log_sum(p for (h, s), p in A.items() if s == s0))
													for s0 in states)
	dic = dict(((s, h), A[(h, s)]-sum_observed[s]) for h, s in A.keys())
	return dic

def update_probabilities(seq, states, hidden_states, pi, Ts, E):
	A = reverse_probs(E, hidden_states, states)
	a, b, xi, gamma = compute_parameters(seq, hidden_states, pi, Ts, A)
	
	E = { }
	for hj in hidden_states:
		cum = log_sum(g[hj] for g in gamma)
		for s in states:
			if sum(1 for s1, g in zip(seq, gamma) if s1 == s) == 0:
				E[(hj,s)] = 0.0
			else:
				E[(hj,s)] = log_sum(g[hj] for s1, g in zip(seq, gamma) if s1 == s) - cum
	
#	total = log_sum(p-E[(h,seq[0])] for h, p in gamma[0].items())
#	pi = dict((h, p-E[(h,seq[0])]-total) for h, p in gamma[0].items())
	pi = gamma[0]
	
	return (pi, E)

def update_pi(seq, states, hidden_states, pi, Ts, E):
	A = reverse_probs(E, hidden_states, states)
	a, b, xi, gamma = compute_parameters(seq, hidden_states, pi, Ts, A)
	pi = gamma[0]
	
	return pi

def print_matrix(matrix, states1, states2):
	def generate_rows():
		yield [''] + states2
		for s1 in states1:
			yield [s1] + [ matrix[(s1,s2)] for s2 in states2 ]
	
	for row in generate_rows():
		print('\t'.join(map(str, row)))

def is_near_parameters(E1, pi1, E2, pi2):
	d2 = log_sum(log_diff(E1[k], E2[k]) for k in E1.keys())
	d3 = log_sum(log_diff(pi1[k], pi2[k]) for k in pi1.keys())
	return log_sum([d2, d3]) < log(1e-7)

def converge_pi_E(Ts, seq, states, hidden_states):
	def is_valid(E):
		return all(v > log(0.5) for k, v in E.items() if k[0] == k[1])
	
	pi = initialize_pi(hidden_states)
	E = initialize_emission_matrix(hidden_states, states)
	for _ in range(100):
		try:
			pi1, E1 = update_probabilities(
								seq, states, hidden_states, pi, Ts, E)
			if is_near_parameters(E, pi, E1, pi1):
				break
			pi, E = pi1, E1
		except ZeroDivisionError:
			break
	
	if is_valid(E):
		return (pi, E)
	else:
		return None

def converge_pi(Ts, seq, states, hidden_states):
	pi = initialize_pi(hidden_states)
	E = initialize_emission_matrix(hidden_states, states)
	for _ in range(100):
		try:
			pi1 = update_pi(seq, states, hidden_states, pi, Ts, E)
			if log_sum(log_diff(pi[k], pi1[k]) for k in pi.keys()) < log(1e-7):
				break
			pi = pi1
		except ZeroDivisionError:
			print(pi)
			break
	
	return (pi, E)

def Baum_Welch(Ts_, seq, states, hidden_states):
	Ts = [ { k: modified_log(v) for k, v in T.items() } for T in Ts_ ]
	result = converge_pi_E(Ts, seq, states, hidden_states)
	if result is not None:
		return result	# (pi, E)
	
	pi, E = converge_pi(Ts, seq, states, hidden_states)
	return (pi, E)


#################### Viterbi ####################

def Viterbi(seq, states, hidden_states, pi, Ts_, A):
	Ts = [ { k: modified_log(v) for k, v in T.items() } for T in Ts_ ]
	
	table = [ { } for _ in seq ]
	table[0] = { h: pi[h]+A[(seq[0],h)] for h in hidden_states }
	
	for k in range(1, len(seq)):
		s = seq[k]
		T = Ts[k-1]
		for h in hidden_states:
			table[k][h] = max(table[k-1][h0]+T[(h0,h)]+A[(s,h)]
											for h0 in hidden_states)
	
	# backtrack
	p, h = max(((table[-1][h], h) for h in hidden_states), key=lambda x: x[0])
	hidden_seq = h
	for k in range(len(seq)-1, 0, -1):
		s = seq[k]
		T = Ts[k-1]
		p, h0 = max((table[k-1][h0]+T[(h0,h)]+A[(s,h)], h0)
											for h0 in hidden_states)
		hidden_seq += h0
		h = h0
	
	return(hidden_seq[::-1])
