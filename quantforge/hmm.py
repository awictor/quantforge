"""Hidden Markov model: forward likelihood and Viterbi decoding.

A discrete HMM has ``n`` hidden states with initial distribution ``pi``, a
state-transition matrix ``A`` (``A[i][j]`` = P(next=j | cur=i)), and emission
matrix ``B`` (``B[i][o]`` = P(observe symbol o | state i)). Given an observation
sequence this module computes:

- ``hmm_forward`` -- the sequence log-likelihood by the scaled forward algorithm
  (scaling each step to avoid underflow), and
- ``hmm_viterbi`` -- the most-likely hidden-state path by dynamic programming in
  log space.

Pure standard library.
"""

import math


def _check(pi, A, B, obs):
    n = len(pi)
    if len(A) != n or any(len(row) != n for row in A):
        raise ValueError("A must be n x n matching pi")
    if len(B) != n:
        raise ValueError("B must have one emission row per state")
    m = len(B[0])
    if any(len(row) != m for row in B):
        raise ValueError("B rows must have equal length")
    if not obs:
        raise ValueError("observation sequence must be non-empty")
    if any(not (0 <= o < m) for o in obs):
        raise ValueError("observation symbol out of range")
    return n


def hmm_forward(pi, A, B, obs):
    """Log-likelihood ``log P(obs | model)`` by the scaled forward algorithm.

    Rescales the forward variables at each step (dividing by their sum) and
    accumulates the log of the scale factors, so the likelihood is exact without
    underflow on long sequences. ``pi`` initial distribution, ``A`` transition, ``B``
    emission matrices; ``obs`` a list of symbol indices.
    """
    n = _check(pi, A, B, obs)
    # Initialize.
    alpha = [pi[i] * B[i][obs[0]] for i in range(n)]
    c = sum(alpha)
    if c <= 0.0:
        return -float("inf")
    alpha = [a / c for a in alpha]
    loglik = math.log(c)
    for t in range(1, len(obs)):
        new = [B[j][obs[t]] * sum(alpha[i] * A[i][j] for i in range(n))
               for j in range(n)]
        c = sum(new)
        if c <= 0.0:
            return -float("inf")
        alpha = [v / c for v in new]
        loglik += math.log(c)
    return loglik


def hmm_viterbi(pi, A, B, obs):
    """Most-likely hidden-state path (Viterbi) and its log-probability.

    Dynamic programming in log space; returns ``(path, log_prob)`` where ``path`` is
    the list of state indices maximizing the joint probability of states and
    observations. Ties are broken toward the lower state index.
    """
    n = _check(pi, A, B, obs)
    T = len(obs)
    neg = -float("inf")

    def ln(x):
        return math.log(x) if x > 0.0 else neg

    # delta[t][j] = best log-prob of a path ending in state j at time t.
    delta = [[neg] * n for _ in range(T)]
    back = [[0] * n for _ in range(T)]
    for j in range(n):
        delta[0][j] = ln(pi[j]) + ln(B[j][obs[0]])
    for t in range(1, T):
        for j in range(n):
            best_i, best_v = 0, neg
            for i in range(n):
                v = delta[t - 1][i] + ln(A[i][j])
                if v > best_v:
                    best_v, best_i = v, i
            delta[t][j] = best_v + ln(B[j][obs[t]])
            back[t][j] = best_i
    # Termination.
    last = max(range(n), key=lambda j: delta[T - 1][j])
    log_prob = delta[T - 1][last]
    path = [last]
    for t in range(T - 1, 0, -1):
        path.append(back[t][path[-1]])
    path.reverse()
    return path, log_prob


def hmm_posterior(pi, A, B, obs):
    """Smoothed posterior state probabilities ``P(state_t = i | obs)`` (gamma).

    Runs the scaled forward-backward algorithm and returns a list of length-``n``
    distributions, one per time step, each summing to one -- the probability of being
    in each hidden state at that time given the *entire* observation sequence (unlike
    Viterbi's single best path, this is the per-time marginal). Pure standard library.
    """
    n = _check(pi, A, B, obs)
    T = len(obs)

    # Scaled forward pass, storing the scale factors c_t.
    alpha = [[0.0] * n for _ in range(T)]
    c = [0.0] * T
    a0 = [pi[i] * B[i][obs[0]] for i in range(n)]
    c[0] = sum(a0) or 1e-300
    alpha[0] = [v / c[0] for v in a0]
    for t in range(1, T):
        row = [B[j][obs[t]] * sum(alpha[t - 1][i] * A[i][j] for i in range(n))
               for j in range(n)]
        c[t] = sum(row) or 1e-300
        alpha[t] = [v / c[t] for v in row]

    # Scaled backward pass with the same c_t.
    beta = [[0.0] * n for _ in range(T)]
    beta[T - 1] = [1.0 / c[T - 1]] * n
    for t in range(T - 2, -1, -1):
        for i in range(n):
            beta[t][i] = sum(A[i][j] * B[j][obs[t + 1]] * beta[t + 1][j]
                             for j in range(n)) / c[t]

    # Posterior gamma_t(i) proportional to alpha_t(i) beta_t(i); normalize per step.
    gamma = []
    for t in range(T):
        g = [alpha[t][i] * beta[t][i] for i in range(n)]
        s = sum(g) or 1e-300
        gamma.append([v / s for v in g])
    return gamma
