##
#  Matchy
#  String searching algorithms for PHP, JavaScript, Python
#
#  @version: 2.0.0
#  https://github.com/foo123/Matchy
#
##
# -*- coding: utf-8 -*-

class Matchy:
    """
    Matchy
    String searching algorithms for PHP, JavaScript, Python
    @version: 1.0.0
    https://github.com/foo123/Matchy
    """

    VERSION = "2.0.0"

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
        suffix = {}
        def is_suffix(k, q, c):
            ks = str(k)
            qs = str(q)
            if ks not in suffix: suffix[ks] = {}
            if qs not in suffix[ks]: suffix[ks][qs] = {}
            if c not in suffix[ks][qs]:
                s1 = pattern[0:k]
                s2 = pattern[0:q] + c
                suffix[ks][qs][c] = (s1 == s2[-len(s1):])
            return suffix[ks][qs][c]

        d = list(map(lambda _: {}, [0] * (m+1)))
        def delta(q, c):
            if c not in in_pattern: return 0
            if c in d[q]: return d[q][c]
            k = min(m, q+1)
            while (0 < k) and not is_suffix(k, q, c): k -= 1
            d[q][c] = k
            return k

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
    def __init__(self, input, type = 'l'):
        if '?' == type: type = {'min':0, 'max':1}
        if '*' == type: type = {'min':0, 'max':INF}
        if '+' == type: type = {'min':1, 'max':INF}
        self.input = input
        self.type = type

    def q0(self):
        type = self.type
        input = self.input
        if isinstance(type, dict):
            if 'min' in type:
                return (input.q0(), 0)
            else:
                return (
                    list(range(0, type['errors']+1, 1)),
                    list(range(0, type['errors']+1, 1))
                )
        if 'l' == type:
            return 0
        if '^' == type:
            return False
        if '$' == type:
            return False
        if '!' == type:
            return input.q0()
        if '|' == type:
            return tuple([nfa.q0() for nfa in input])
        if ',' == type:
            return (input[0].q0(), 0)

    def d(self, q, c):
        type = self.type
        input = self.input
        if isinstance(type, dict):
            if 'min' in type:
                q = (input.d(q[0], c), q[1])
                if input.accept(q[0]): q = (input.q0(), q[1]+1)
                return q
            else:
                k = type['errors']
                w = input
                n = len(w)
                index, value = q
                new_index = []
                new_value = []
                m = len(index)
                last_i = -1
                last_v = 0
                next_i = -1
                if (0 < m) and (0 == index[0]) and (value[0] < k):
                    last_i = 0
                    last_v = value[0] + 1
                    new_index.append(last_i)
                    new_value.append(last_v)

                for j, i in enumerate(index):
                    if i >= n: break
                    d = 0 if w[i] == c else 1
                    v = value[j] + d
                    next_i = index[j+1] if j+1 < m else -1
                    if i == last_i:
                        v = min(v, last_v + 1)
                    if i+1 == next_i:
                        v = min(v, value[j+1] + 1)
                    if v <= k:
                        last_i = i+1
                        last_v = v
                        new_index.append(last_i)
                        new_value.append(last_v)
                return (
                    new_index,
                    new_value
                )
        if 'l' == type:
            if isinstance(c, int): return q
            return (q+1 if c == input[q] else 0) if q < len(input) else 0
        if '^' == type:
            return 0 == c
        if '$' == type:
            return 1 == c
        if '!' == type:
            return input.d(q, c);
        if '|' == type:
            return tuple([nfa.d(q[i], c) for i, nfa in enumerate(input)])
        if ',' == type:
            i = q[1]
            if input[i].accept(q[0]):
                if i+1 < len(input):
                    q0 = (input[i].d(q[0], c), i)
                    q1 = (input[i+1].d(input[i+1].q0(), c), i+1)
                    return q0 if input[i+1].reject(q1[0]) and not input[i].reject(q0[0]) else q1
                else:
                    return q
            else:
                return (input[i].d(q[0], c), i)

    def accept(self, q):
        type = self.type
        input = self.input
        if isinstance(type, dict):
            if 'min' in type:
                return (type['min'] <= q[1]) and (q[1] <= type['max'])
            else:
                return (0 < len(q[0])) and (q[0][-1] >= len(input))
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

    def reject(self, q):
        type = self.type
        input = self.input
        if isinstance(type, dict):
            if 'min' in type:
                return ((q[1] >= type['max']) and input.accept(q[0])) or ((q[1] < type['min']) and input.reject(q[0]))
            else:
                return not q[0]
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

    def match(self, string, offset = 0, q = None):
        i = offset
        j = i
        n = len(string)
        if q is None: q = self.q0()
        while True:
            if j >= n:
                i += 1
                if i >= n: break
                j = i
                q = self.q0()
            if 0 == j: q = self.d(q, 0)
            q = self.d(q, string[j])
            if n-1 == j: q = self.d(q, 1)
            if self.accept(q):
                return i # matched
            elif self.reject(q):
                j = n # failed, try next
            else:
                j += 1 # continue
        return -1

Matchy.NFA = NFA

INF = float('inf')

__all__ = ['Matchy']

