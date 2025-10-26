##
#  Matchy
#  Exact and fuzzy string searching algorithms for PHP, JavaScript, Python
#
#  @version: 3.0.0
#  https://github.com/foo123/Matchy
#
##
# -*- coding: utf-8 -*-
import math

class Matchy:
    """
    Matchy
    String searching algorithms for PHP, JavaScript, Python
    @version: 1.0.0
    https://github.com/foo123/Matchy
    """

    VERSION = "3.0.0"

    def __init__(self):
        pass

    def fsa(self, pattern, string = None, offset = 0):
        # https://en.wikipedia.org/wiki/Finite-state_machine
        # https://en.wikipedia.org/wiki/String-searching_algorithm
        # https://euroinformatica.ro/documentation/programming/!!!Algorithms_CORMEN!!!/DDU0214.html

        # construct transition matrix
        m = len(pattern)
        in_pattern = {}
        for i in range(m):
            in_pattern[pattern[i]] = 1
        def is_suffix(k, q, c):
            s1 = pattern[0:k]
            s2 = pattern[0:q] + c
            return (s1 == s2[-len(s1):])

        d = list(map(lambda _: {}, [0] * (m+1)))
        def delta(q, c):
            if c not in in_pattern: return 0
            if c not in d[q]:
                k = min(m, q+1)
                while (0 < k) and not is_suffix(k, q, c): k -= 1
                d[q][c] = k
            return d[q][c]

        # matcher
        def matcher(string, offset = 0):
            n = len(string)
            if 0 > offset: offset += n
            if (0 < n) and (0 < m) and (n >= offset+m):
                q = 0
                for i in range(offset, n):
                    c = string[i]
                    q = delta(q, c)
                    if m == q: return i-m+1 # matched
            return -1

        return matcher if string is None else matcher(string, offset)

    def rabinkarp(self, pattern, string = None, offset = 0):
        # https://en.wikipedia.org/wiki/Rabin%E2%80%93Karp_algorithm
        # http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=BA184C94C16CB23D5FA7329E257E3713?doi=10.1.1.86.9502&rep=rep1&type=pdf

        m = len(pattern)
        Q = 3989 # prime number so that 10*Q can fit in one WORD (i.e 2^16 bits)

        # matcher
        def matcher(string, offset = 0):
            n = len(string)
            if 0 > offset: offset += n
            if (0 < n) and (0 < m) and (n >= offset+m):
                alphabet = {}
                D = 0
                for i in range(m):
                    c = pattern[i]
                    if c not in alphabet:
                        alphabet[c] = D
                        D += 1
                for i in range(offset, n):
                    c = string[i]
                    if c not in alphabet:
                        alphabet[c] = D
                        D += 1
                h = 1
                pq = alphabet[pattern[0]] % Q
                sq = alphabet[string[offset]] % Q
                for i in range(1, m):
                    h  = (h*D) % Q
                    pq = (pq*D + alphabet[pattern[i]]) % Q
                    sq = (sq*D + alphabet[string[offset+i]]) % Q

                n = n-m
                for i in range(offset, n+1):
                    # pq, sq, D-base arithmetic code, used as a quick "hash" test
                    # worst case: many hash collisions -> "naive" matching
                    if pq == sq:
                        if string[i:i+m] == pattern: return i # matched
                    # update text hash for next char using Horner algorithm
                    if i < n:
                        sq = (D*(sq - h*alphabet[string[i]]) + alphabet[string[i+m]]) % Q
                        if 0 > sq: sq += Q
            return -1

        return matcher if string is None else matcher(string, offset)

    def knuthmorrispratt(self, pattern, string = None, offset = 0):
        # https://en.wikipedia.org/wiki/Knuth%E2%80%93Morris%E2%80%93Pratt_algorithm
        # http://www.eecs.ucf.edu/~shzhang/Combio09/kmp.pdf

        # pre-processing
        m = len(pattern)
        T = [0] * m
        T[0] = -1
        k = 1
        i = 0
        while k < m:
            c = pattern[k]
            if c == pattern[i]:
                T[k] = T[i]
            else:
                T[k] = i
                while (0 <= i) and (c != pattern[i]): i = T[i]
            k += 1
            i += 1

        # matcher
        def matcher(string, offset = 0):
            n = len(string)
            if 0 > offset: offset += n
            if (0 < n) and (0 < m) and (n >= offset+m):
                k = 0
                i = offset
                while i < n:
                    if pattern[k] == string[i]:
                        i += 1
                        k += 1
                        if k == m: return i - k # matched
                    else:
                        k = T[k]
                        if k < 0:
                            i += 1
                            k += 1
            return -1

        return matcher if string is None else matcher(string, offset)

    def boyermoore(self, pattern, string = None, offset = 0):
        # https://en.wikipedia.org/wiki/Boyer%E2%80%93Moore_string-search_algorithm
        # http://www.cs.utexas.edu/~moore/publications/fstrpos.pdf

        # pre-processing
        m = len(pattern)
        suffix = [0] * (m+1)
        shift = [0] * (m+1)
        i = m
        j = m+1
        suffix[i] = j
        while 0 < i:
            while (j <= m) and (pattern[i-1] != pattern[j-1]):
                if 0 == shift[j]: shift[j] = j-i
                j = suffix[j]
            i -= 1
            j -= 1
            suffix[i] = j
        j = suffix[0]
        for i in range(m):
            if 0 == shift[i]:
                shift[i] = j
                if i == j: j = suffix[j]

        # matcher
        def matcher(string, offset = 0):
            n = len(string)
            if 0 > offset: offset += n
            if (0 < n) and (0 < m) and (n >= offset+m):
                i = offset
                while i + m <= n:
                    j = m - 1
                    while (0 <= j) and (pattern[j] == string[i+j]):
                        j -= 1
                    if 0 > j:
                        return i # matched
                    else:
                        i += shift[j+1]
            return -1

        return matcher if string is None else matcher(string, offset)

    def twoway(self, pattern, string = None, offset = 0):
        # https://en.wikipedia.org/wiki/Two-way_string-matching_algorithm
        # http://monge.univ-mlv.fr/~mac/Articles-PDF/CP-1991-jacm.pdf

        # pre-processing
        m = len(pattern)
        def maximal_suffix(s, n, dir):
            # assumes 1-indexed strings instead of 0-indexed
            p = 1  #  currently known period.
            k = 1  # index for period testing, 0 < k <= p.
            j = 1  # index for maxsuf testing. greater than maxs.
            i = 0  # the proposed starting index of maxsuf

            while j + k <= n:
                a = s[j + k - 1]
                b = s[i + k - 1]
                cmp = dir*(1 if a > b else (-1 if a < b else 0))
                if 0 > cmp:
                    # Suffix is smaller. Period is the entire prefix so far.
                    j += k
                    k = 1
                    p = j - i
                elif 0 < cmp:
                    # Suffix is larger. Start over from here.
                    i = j
                    j = i + 1
                    k = 1
                    p = 1
                else:
                    # They are the same - we should go on.
                    if k == p:
                        # We are done checking this stretch of p. reset k.
                        j += p
                        k = 1
                    else:
                        k += 1

            return (i, p)

        # critical factorization
        s1 = maximal_suffix(pattern, m,  1)
        s2 = maximal_suffix(pattern, m, -1)
        if s1[0] >= s2[0]:
            l, p = s1
        else:
            l, p = s2
        # small period
        small_period = (2*l < m) and (pattern[0:l] == pattern[l:l+p][-l:])

        # matcher
        def matcher(string, offset = 0):
            n = len(string)
            if 0 > offset: offset += n
            if (0 < n) and (0 < m) and (n >= offset+m):
                pos = offset
                if small_period:
                    # small period
                    s = 0
                    while pos + m <= n:
                        i = max(l, s) + 1
                        while (i <= m) and (pattern[i-1] == string[pos+i-1]):
                            i += 1
                        if i <= m:
                            pos += max(i-l, s-p+1)
                            s = 0
                        else:
                            j = l
                            while (j > s) and (pattern[j-1] == string[pos+j-1]):
                                j -= 1
                            if j <= s: return pos # matched
                            pos += p
                            s = m-p
                else:
                    # large period
                    q = max(l, m-l) + 1
                    while pos + m <= n:
                        i = l + 1
                        while (i <= m) and (pattern[i-1] == string[pos+i-1]):
                            i += 1
                        if i <= m:
                            pos += i-l
                        else:
                            j = l
                            while (j > 0) and (pattern[j-1] == string[pos+j-1]):
                                j -= 1
                            if 0 == j: return pos # matched
                            pos += q
            return -1

        return matcher if string is None else matcher(string, offset)

    def commentzwalter(self, pattern, string = None, offset = 0):
        # https://en.wikipedia.org/wiki/Commentz-Walter_algorithm
        return non_matcher if string is None else non_matcher(string, offset) # TODO

    def baezayatesgonnet(self, pattern, string = None, offset = 0):
        # https://en.wikipedia.org/wiki/Bitap_algorithm
        return non_matcher if string is None else non_matcher(string, offset) # TODO

    def ahocorasick(self, pattern, string = None, offset = 0):
        # https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_algorithm
        return non_matcher if string is None else non_matcher(string, offset) # TODO


def non_matcher(string, offset = 0):
    return -1

# non-deterministic finite automaton
class NFA:
    def makeFuzzy(nfa, max_errors, transpositions = False):
        if isinstance(nfa, str):
            # change literal string to approximate match
            nfa = NFA(nfa, {'errors':max_errors, 'transpositions':transpositions})
        elif isinstance(nfa, NFA):
            if isinstance(nfa.type, dict) and ('errors' in nfa.type):
                # leave as is
                pass
            elif 'l' == nfa.type:
                # change literal match to approximate match
                nfa = NFA(nfa.input, {'errors':max_errors, 'transpositions':transpositions})
            elif isinstance(nfa.input, (list,tuple)):
                # recurse
                nfa.input = list(map(lambda input: NFA.makeFuzzy(input, max_errors, transpositions), nfa.input))
            elif isinstance(nfa.input, NFA):
                # recurse
                nfa.input = NFA.makeFuzzy(nfa.input, max_errors, transpositions)
        return nfa

    def __init__(self, input, type = 'l'):
        if '?' == type: type = {'min':0, 'max':1}
        if '*' == type: type = {'min':0, 'max':INF}
        if '+' == type: type = {'min':1, 'max':INF}
        self.type = type
        self.input = NFA.makeFuzzy(input, type['total_errors'], type['transpositions'] if 'transpositions' in type else False) if isinstance(type, dict) and ('total_errors' in type) else input

    def q0(self):
        type = self.type
        input = self.input
        if isinstance(type, dict):
            if 'min' in type:
                q = (input.q0(), 0)
            elif 'errors' in type:
                k = min(type['errors'], len(input))
                q = (
                    list(range(0, k+1, 1)),
                    list(range(0, k+1, 1)),
                    [],
                    [],
                    ''
                ) if ('transpositions' in type) and (type['transpositions']) else (
                    list(range(0, k+1, 1)),
                    list(range(0, k+1, 1))
                )
            else:
                q = input.q0()
        if 'l' == type:
            q = 0
        if '^' == type:
            q = False
        if '$' == type:
            q = False
        if '!' == type:
            q = input.q0()
        if '|' == type:
            q = tuple([nfa.q0() for nfa in input])
        if ',' == type:
            q = (input[0].q0(), 0)
        return {'q':q, 'e':0} # keep track of errors

    def d(self, qe, c):
        type = self.type
        input = self.input
        q = qe['q']
        e = qe['e']
        if isinstance(type, dict):
            if 'min' in type:
                e0 = q[0]['e']
                qe = input.d(q[0], c)
                q = (qe, q[1])
                e += q[0]['e'] - e0
                if input.accept(qe): q = (input.q0(), q[1]+1)
            elif 'errors' in type:
                if isinstance(c, int):
                    #q = q
                    pass
                else:
                    transpositions = ('transpositions' in type) and (type['transpositions'])
                    w = input
                    n = len(w)
                    k = min(type['errors'], n)
                    min_e = k+1
                    index = q[0]
                    value = q[1]
                    new_index = []
                    new_value = []
                    m = len(index)
                    prev_i = -1
                    prev_v = 0
                    next_i = -1
                    if transpositions:
                        index_2 = q[2]
                        value_2 = q[3]
                        cp = q[4]
                        m2 = len(index_2)
                        j2 = 0
                    if (0 < m) and (0 == index[0]) and (value[0] < k):
                        i = 0
                        v = value[0] + 1
                        prev_i = i
                        prev_v = v
                        new_index.append(i)
                        new_value.append(v)

                    for j, i in enumerate(index):
                        if i >= n: break
                        d = 0 if w[i] == c else 1
                        v = value[j] + d # L[i,ii] = L[i-1,ii-1] + d
                        next_i = index[j+1] if j+1 < m else -1
                        i += 1
                        if i-1 == prev_i:
                            v = min(v, prev_v + 1) # L[i,ii] = min(L[i,ii], L[i-1,ii] + 1)
                        if i == next_i:
                            v = min(v, value[j+1] + 1) # L[i,ii] = min(L[i,ii], L[i,ii-1] + 1)
                        if transpositions and (cp == w[i-1]) and (c == w[i-2]):
                            while (j2 < m2) and (index_2[j2] < i-2): j2 += 1
                            if (j2 < m2) and (i-2 == index_2[j2]):
                                v = min(v, value_2[j2] + d) # L[i,ii] = min(L[i,ii], L[i-2,ii-2] + d)
                                j2 += 1
                        if v <= k:
                            min_e = min(min_e, v)
                            prev_i = i
                            prev_v = v
                            new_index.append(i)
                            new_value.append(v)
                    q = (
                        new_index,
                        new_value,
                        index,
                        value,
                        c
                    ) if transpositions else (
                        new_index,
                        new_value
                    )
                    e = min_e
            else:
                q = input.d(q, c)
                e = q['e']
        if 'l' == type:
            if isinstance(c, int):
                # q = q
                pass
            else:
                q = (q+1 if c == input[q] else 0) if q < len(input) else 0
        if '^' == type:
            q = isinstance(c, int) and (0 == c)
            #e = 0 if q else 1
        if '$' == type:
            q = isinstance(c, int) and (1 == c)
            #e = 0 if q else 1
        if '!' == type:
            q = input.d(q, c);
            e = q['e']
        if '|' == type:
            e0 = e
            e = INF
            q = tuple([nfa.d(q[i], c) for i, nfa in enumerate(input)])
            for i, qi in enumerate(q):
                if not input[i].reject(qi): e = min(e, qi['e'])
            if not math.isfinite(e): e = e0+1
        if ',' == type:
            i = q[1]
            if input[i].accept(q[0]):
                if i+1 < len(input):
                    e0 = q[0]['e']
                    q0 = input[i].d(q[0], c)
                    q1 = input[i+1].d(input[i+1].q0(), c)
                    if input[i+1].reject(q1) and not input[i].reject(q0):
                        q = (q0, i)
                        e += q0['e'] - e0
                    else:
                        q = (q1, i+1)
                        e += q1['e']
                else:
                    #q = q
                    pass
            else:
                e0 = q[0]['e']
                q = (input[i].d(q[0], c), i)
                e += q[0]['e'] - e0
        return {'q':q, 'e':e} # keep track of errors

    def accept(self, qe):
        type = self.type
        input = self.input
        q = qe['q']
        e = qe['e']
        if isinstance(type, dict):
            if 'min' in type:
                return (type['min'] <= q[1]) and (q[1] <= type['max'])
            elif 'errors' in type:
                return (0 < len(q[0])) and (q[0][-1] == len(input))
            else:
                return (e <= type['total_errors']) and input.accept(q)
        if 'l' == type:
            return q == len(input)
        if '^' == type:
            return q
        if '$' == type:
            return q
        if '!' == type:
            return input.reject(q)
        if '|' == type:
            return 0 < len(list(filter(lambda entry: entry[1].accept(q[entry[0]]), enumerate(input))))
        if ',' == type:
            return (q[1]+1 == len(input)) and input[q[1]].accept(q[0])

    def reject(self, qe):
        type = self.type
        input = self.input
        q = qe['q']
        e = qe['e']
        if isinstance(type, dict):
            if 'min' in type:
                return ((q[1] >= type['max']) and input.accept(q[0])) or ((q[1] < type['min']) and input.reject(q[0]))
            elif 'errors' in type:
                return not q[0]
            else:
                return (e > type['total_errors']) or input.reject(q)
        if 'l' == type:
            return 0 == q
        if '^' == type:
            return not q
        if '$' == type:
            return not q
        if '!' == type:
            return input.accept(q)
        if '|' == type:
            return len(input) == len(list(filter(lambda entry: entry[1].reject(q[entry[0]]), enumerate(input))))
        if ',' == type:
            return input[q[1]].reject(q[0])

    def match(self, string, offset = 0, return_match = False, q = None):
        i = offset
        j = i
        n = len(string)
        if q is None: q = self.q0()
        c = ''
        while True:
            if j >= n:
                i += 1
                if i >= n: break
                j = i
                q = self.q0()
                c = string[j-1]
            prevc = c
            c = string[j]
            if (0 == j) or ("\n" == prevc): q = self.d(q, 0)
            q = self.d(q, c)
            if (n-1 == j) or ("\n" == c): q = self.d(q, 1)
            if self.accept(q):
                return string[i:j+1] if return_match else i # matched
            elif self.reject(q):
                j = n # failed, try next
            else:
                j += 1 # continue
        return None if return_match else -1

Matchy.NFA = NFA

INF = float('inf')

__all__ = ['Matchy']

