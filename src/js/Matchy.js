/**
*  Matchy
*  Exact and fuzzy string searching algorithms for PHP, JavaScript, Python
*
*  @version: 2.0.0
*  https://github.com/foo123/Matchy
*
**/
!function(root, name, factory) {
"use strict";
if (('object' === typeof module) && module.exports) /* CommonJS */
    (module.$deps = module.$deps||{}) && (module.exports = module.$deps[name] = factory.call(root));
else if (('function' === typeof define) && define.amd && ('function' === typeof require) && ('function' === typeof require.specified) && require.specified(name) /*&& !require.defined(name)*/) /* AMD */
    define(name, ['module'], function(module) {factory.moduleUri = module.uri; return factory.call(root);});
else if (!(name in root)) /* Browser/WebWorker/.. */
    (root[name] = factory.call(root)||1) && ('function' === typeof(define)) && define.amd && define(function() {return root[name];});
}(/* current root */          'undefined' !== typeof self ? self : this,
  /* module name */           "Matchy",
  /* module factory */        function ModuleFactory__Matchy(undef) {
"use strict";

function Matchy()
{
}
Matchy.VERSION = "2.0.0";
Matchy.prototype = {
    constructor: Matchy,

    fsa: function(pattern, string, offset) {
        // https://en.wikipedia.org/wiki/Finite-state_machine
        // https://en.wikipedia.org/wiki/String-searching_algorithm
        // https://euroinformatica.ro/documentation/programming/!!!Algorithms_CORMEN!!!/DDU0214.html

        // construct transition matrix
        var m = pattern.length;
        var delta, d, is_suffix, in_pattern;

        in_pattern = {};
        for (var i=0; i<m; ++i)
        {
            in_pattern[pattern.charAt(i)] = 1;
        }
        is_suffix = function(k, q, c) {
            var s1 = pattern.slice(0, k),
                s2 = pattern.slice(0, q) + c;
            return (s1 === s2.slice(-s1.length));
        };

        d = array_fill(0, m+1, 0).map(function() {return {};});
        delta = function(q, c) {
            if (!isset(in_pattern, c) || (1 !== in_pattern[c])) return 0;
            if (isset(d[q], c)) return d[q][c];
            var k = min(m, q+1);
            while ((0 < k) && !is_suffix(k, q, c)) --k;
            d[q][c] = k;
            return k;
        };

        // matcher
        var matcher = function(string, offset) {
            var n = string.length;
            if (null == offset) offset = 0;
            if (0 > offset) offset += n;
            if ((0 < n) && (0 < m) && (n >= offset+m))
            {
                for (var c,q=0,i=offset; i<n; ++i)
                {
                    c = string.charAt(i);
                    q = delta(q, c);
                    if (m === q) return i-m+1; // matched
                }
            }
            return -1;
        };
        return null == string ? matcher : matcher(string, offset);
    },

    rabinkarp: function(pattern, string, offset) {
        // https://en.wikipedia.org/wiki/Rabin%E2%80%93Karp_algorithm
        // http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=BA184C94C16CB23D5FA7329E257E3713?doi=10.1.1.86.9502&rep=rep1&type=pdf

        var m = pattern.length;
        var Q = 3989; // prime number so that 10*Q can fit in one WORD (i.e 2^16 bits)

        // matcher
        var matcher = function(string, offset) {
            var n = string.length;
            if (null == offset) offset = 0;
            if (0 > offset) offset += n;
            if ((0 < n) && (0 < m) && (n >= offset+m))
            {
                var h, pq, sq, i, c,
                    alphabet = {}, D = 0
                ;
                for (i=0; i<m; ++i)
                {
                    c = pattern.charAt(i);
                    if (!isset(alphabet, c)) alphabet[c] = D++;
                }
                for (i=offset; i<n; ++i)
                {
                    c = string.charAt(i);
                    if (!isset(alphabet, c)) alphabet[c] = D++;
                }
                h = 1;
                pq = alphabet[pattern.charAt(0)] % Q;
                sq = alphabet[string.charAt(offset)] % Q;
                for (i=1; i<m; ++i)
                {
                    h  = (h*D) % Q;
                    pq = (pq*D + alphabet[pattern.charAt(i)]) % Q;
                    sq = (sq*D + alphabet[string.charAt(offset+i)]) % Q;
                }

                n = n-m;
                for (i=offset; i<=n; ++i)
                {
                    // pq, sq, D-base arithmetic code, used as a quick "hash" test
                    // worst case: many hash collisions -> "naive" matching
                    if (pq === sq)
                    {
                        if (string.slice(i, i+m) === pattern) return i; // matched
                    }
                    // update text hash for next char using Horner algorithm
                    if (i < n)
                    {
                        sq = (D*(sq - h*alphabet[string.charAt(i)]) + alphabet[string.charAt(i+m)]) % Q;
                        if (0 > sq) sq += Q;
                    }
                }
            }
            return -1;
        };
        return null == string ? matcher : matcher(string, offset);
    },

    knuthmorrispratt: function(pattern, string, offset) {
        // https://en.wikipedia.org/wiki/Knuth%E2%80%93Morris%E2%80%93Pratt_algorithm
        // http://www.eecs.ucf.edu/~shzhang/Combio09/kmp.pdf

        // pre-processing
        var m = pattern.length;
        var i, k, c, T = array_fill(0, m, 0);
        T[0] = -1;
        k = 1;
        i = 0;
        while (k < m)
        {
            c = pattern.charAt(k);
            if (c === pattern.charAt(i))
            {
                T[k] = T[i];
            }
            else
            {
                T[k] = i;
                while ((0 <= i) && (c !== pattern.charAt(i))) i = T[i];
            }
            ++k;
            ++i;
        }

        // matcher
        var matcher = function(string, offset) {
            var n = string.length;
            if (null == offset) offset = 0;
            if (0 > offset) offset += n;
            if ((0 < n) && (0 < m) && (n >= offset+m))
            {
                var k = 0, i = offset;
                while (i < n)
                {
                    if (pattern.charAt(k) === string.charAt(i))
                    {
                        ++i;
                        ++k;
                        if (k === m) return i - k; // matched
                    }
                    else
                    {
                        k = T[k];
                        if (k < 0)
                        {
                            ++i;
                            ++k;
                        }
                    }
                }
            }
            return -1;
        };
        return null == string ? matcher : matcher(string, offset);
    },

    boyermoore: function(pattern, string, offset) {
        // https://en.wikipedia.org/wiki/Boyer%E2%80%93Moore_string-search_algorithm
        // http://www.cs.utexas.edu/~moore/publications/fstrpos.pdf

        // pre-processing
        var m = pattern.length;
        var suffix = array_fill(0, m+1, 0),
            shift = array_fill(0, m+1, 0),
            i, j;
        i = m;
        j = m+1;
        suffix[i] = j;
        while (0 < i)
        {
            while ((j <= m) && (pattern.charAt(i-1) !== pattern.charAt(j-1)))
            {
                if (0 === shift[j]) shift[j] = j-i;
                j = suffix[j];
            }
            --i;
            --j;
            suffix[i] = j;
        }
        j = suffix[0];
        for (i=0; i<m; ++i)
        {
            if (0 === shift[i])
            {
                shift[i] = j;
                if (i === j) j = suffix[j];
            }
        }

        // matcher
        var matcher = function(string, offset) {
            var n = string.length;
            if (null == offset) offset = 0;
            if (0 > offset) offset += n;
            if ((0 < n) && (0 < m) && (n >= offset+m))
            {
                var i = offset, j;
                while (i + m <= n)
                {
                    j = m - 1;
                    while ((0 <= j) && (pattern.charAt(j) === string.charAt(i+j)))
                    {
                        --j;
                    }
                    if (0 > j)
                    {
                        return i; // matched
                    }
                    else
                    {
                        i += shift[j+1];
                    }
                }
            }
            return -1;
        };
        return null == string ? matcher : matcher(string, offset);
    },

    twoway: function(pattern, string, offset) {
        // https://en.wikipedia.org/wiki/Two-way_string-matching_algorithm
        // http://monge.univ-mlv.fr/~mac/Articles-PDF/CP-1991-jacm.pdf

        // pre-processing
        var m = pattern.length;
        var maximal_suffix = function(s, n, dir) {
            // assumes 1-indexed strings instead of 0-indexed
            var p = 1,  //  currently known period.
                k = 1,  // index for period testing, 0 < k <= p.
                j = 1,  // index for maxsuf testing. greater than maxs.
                i = 0,  // the proposed starting index of maxsuf
                a, b, cmp
            ;
            while (j + k <= n)
            {
                a = s.charAt(j + k - 1);
                b = s.charAt(i + k - 1);
                cmp = dir*(a > b ? 1 : (a < b ? -1 : 0));
                if (0 > cmp)
                {
                    // Suffix is smaller. Period is the entire prefix so far.
                    j += k;
                    k = 1;
                    p = j - i;
                }
                else if (0 < cmp)
                {
                    // Suffix is larger. Start over from here.
                    i = j;
                    j = i + 1;
                    k = 1;
                    p = 1;
                }
                else
                {
                    // They are the same - we should go on.
                    if (k === p)
                    {
                        // We are done checking this stretch of p. reset k.
                        j += p;
                        k = 1;
                    }
                    else
                    {
                        ++k;
                    }
                }
            }
            return [i, p];
        };
        // critical factorization
        var s1, s2, l, p;
        s1 = maximal_suffix(pattern, m,  1);
        s2 = maximal_suffix(pattern, m, -1);
        if (s1[0] >= s2[0])
        {
            l = s1[0];
            p = s1[1];
        }
        else
        {
            l = s2[0];
            p = s2[1];
        }
        // small period
        var small_period = (2*l < m) && (pattern.slice(0, l) === pattern.slice(l, l+p).slice(-l));

        // matcher
        var matcher = function(string, offset) {
            var n = string.length;
            if (null == offset) offset = 0;
            if (0 > offset) offset += n;
            if ((0 < n) && (0 < m) && (n >= offset+m))
            {
                // assumes 1-indexed strings instead of 0-indexed
                var pos, i, j, s, q;

                pos = offset;
                if (small_period)
                {
                    // small period
                    s = 0;
                    while (pos + m <= n)
                    {
                        i = max(l, s) + 1;
                        while ((i <= m) && (pattern.charAt(i-1) === string.charAt(pos+i-1)))
                        {
                            ++i;
                        }
                        if (i <= m)
                        {
                            pos += max(i-l, s-p+1);
                            s = 0;
                        }
                        else
                        {
                            j = l;
                            while ((j > s) && (pattern.charAt(j-1) === string.charAt(pos+j-1)))
                            {
                                --j;
                            }
                            if (j <= s) return pos; // matched
                            pos += p;
                            s = m-p;
                        }
                    }
                }
                else
                {
                    // large period
                    q = max(l, m-l) + 1;
                    while (pos + m <= n)
                    {
                        i = l + 1;
                        while ((i <= m) && (pattern.charAt(i-1) === string.charAt(pos+i-1)))
                        {
                            ++i;
                        }
                        if (i <= m)
                        {
                            pos += i-l;
                        }
                        else
                        {
                            j = l;
                            while ((j > 0) && (pattern.charAt(j-1) === string.charAt(pos+j-1)))
                            {
                                --j;
                            }
                            if (0 === j) return pos; // matched
                            pos += q;
                        }
                    }
                }
            }
            return -1;
        };
        return null == string ? matcher : matcher(string, offset);
    },

    commentzwalter: function(pattern, string, offset) {
        // https://en.wikipedia.org/wiki/Commentz-Walter_algorithm
        return null == string ? non_matcher : non_matcher(string, offset); // TODO
    },

    baezayatesgonnet: function(pattern, string, offset) {
        // https://en.wikipedia.org/wiki/Bitap_algorithm
        return null == string ? non_matcher : non_matcher(string, offset); // TODO
    },

    ahocorasick: function(pattern, string, offset) {
        // https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_algorithm
        return null == string ? non_matcher : non_matcher(string, offset); // TODO
    }
};

function non_matcher(string, offset)
{
    return -1;
}

// non-deterministic finite automaton
function NFA(input, type)
{
    var self = this;
    if (!(self instanceof NFA)) return new NFA(input, type);
    if (null == type) type = 'l';
    if ('?' === type) type = {min:0, max:1};
    if ('*' === type) type = {min:0, max:INF};
    if ('+' === type) type = {min:1, max:INF};
    self.input = input;
    self.type = type;
}
NFA.prototype = {
    constructor: NFA,

    input: null,
    type: null,

    q0: function() {
        var self = this, type = self.type, input = self.input;
        if (is_obj(type))
        {
            if (isset(type, 'min'))
            {
                return [input.q0(), 0];
            }
            else
            {
                var k = min(type.errors, input.length);
                return type.transpositions ? [
                    range(0, k, 1),
                    range(0, k, 1),
                    [],
                    [],
                    ''
                ] : [
                    range(0, k, 1),
                    range(0, k, 1)
                ];
            }
        }
        if ('l' === type)
        {
            return 0;
        }
        if ('^' === type)
        {
            return false;
        }
        if ('$' === type)
        {
            return false;
        }
        if ('!' === type)
        {
            return input.q0();
        }
        if ('|' === type)
        {
            return input.map(function(nfa) {return nfa.q0();});
        }
        if (',' === type)
        {
            return [input[0].q0(), 0];
        }
    },

    d: function(q, c) {
        var self = this, type = self.type, input = self.input;
        if (is_obj(type))
        {
            if (isset(type, 'min'))
            {
                q = [input.d(q[0], c), q[1]];
                if (input.accept(q[0])) q = [input.q0(), q[1]+1];
                return q;
            }
            else
            {
                if (is_number(c)) return q;
                var transpositions = !!type.transpositions,
                    w = input,
                    n = input.length,
                    k = min(type.errors, n),
                    index = q[0],
                    value = q[1],
                    index_2 = null,
                    value_2 = null,
                    new_index = [],
                    new_value = [],
                    m = index.length,
                    m2 = 0,
                    prev_i = -1,
                    prev_v = 0,
                    next_i = -1,
                    i, j, j2,
                    v, d, cp
                ;
                if (transpositions)
                {
                    index_2 = q[2];
                    value_2 = q[3];
                    cp = q[4];
                    m2 = index_2.length;
                    j2 = 0;
                }
                if ((0 < m) && (0 === index[0]) && (value[0] < k))
                {
                    i = 0;
                    v = value[0] + 1;
                    prev_i = i;
                    prev_v = v;
                    new_index.push(i);
                    new_value.push(v);
                }
                for (j=0; j<m; ++j)
                {
                    i = index[j];
                    if (i >= n) break;
                    d = w.charAt(i) === c ? 0 : 1;
                    v = value[j] + d; // L[i,ii] = L[i-1,ii-1] + d
                    next_i = j+1 < m ? index[j+1] : -1;
                    ++i;
                    if (i-1 === prev_i)
                    {
                        v = min(v, prev_v + 1); // L[i,ii] = min(L[i,ii], L[i-1,ii] + 1)
                    }
                    if (i === next_i)
                    {
                        v = min(v, value[j+1] + 1); // L[i,ii] = min(L[i,ii], L[i,ii-1] + 1)
                    }
                    if (transpositions && (cp === w.charAt(i-1)) && (c === w.charAt(i-2)))
                    {
                        while ((j2 < m2) && (index_2[j2] < i-2)) ++j2;
                        if ((j2 < m2) && (i-2 === index_2[j2]))
                        {
                            v = min(v, value_2[j2] + d); // L[i,ii] = min(L[i,ii], L[i-2,ii-2] + d)
                        }
                    }
                    if (v <= k)
                    {
                        prev_i = i;
                        prev_v = v;
                        new_index.push(i);
                        new_value.push(v);
                    }
                }
                return transpositions ? [
                    new_index,
                    new_value,
                    index,
                    value,
                    c
                ] : [
                    new_index,
                    new_value
                ];
            }
        }
        if ('l' === type)
        {
            if (is_number(c)) return q;
            return q < input.length ? (c === input.charAt(q) ? q+1 : 0) : 0;
        }
        if ('^' === type)
        {
            return is_number(c) && (0 === c);
        }
        if ('$' === type)
        {
            return is_number(c) && (1 === c);
        }
        if ('!' === type)
        {
            return input.d(q, c);
        }
        if ('|' === type)
        {
            return input.map(function(nfa, i) {return nfa.d(q[i], c);});
        }
        if (',' === type)
        {
            var i = q[1], q0, q1;
            if (input[i].accept(q[0]))
            {
                if (i+1 < input.length)
                {
                    q0 = [input[i].d(q[0], c), i];
                    q1 = [input[i+1].d(input[i+1].q0(), c), i+1];
                    return input[i+1].reject(q1[0]) && !input[i].reject(q0[0]) ? q0 : q1;
                }
                else
                {
                    return q;
                }
            }
            else
            {
                return [input[i].d(q[0], c), i];
            }
        }
    },

    accept: function(q) {
        var self = this, type = self.type, input = self.input;
        if (is_obj(type))
        {
            if (isset(type, 'min'))
            {
                return (type.min <= q[1]) && (q[1] <= type.max);
            }
            else
            {
                return (0 < q[0].length) && (q[0][q[0].length-1] === input.length);
            }
        }
        if ('l' === type)
        {
            return q === input.length;
        }
        if ('^' === type)
        {
            return q;
        }
        if ('$' === type)
        {
            return q;
        }
        if ('!' === type)
        {
            return input.reject(q);
        }
        if ('|' === type)
        {
            return 0 < input.filter(function(nfa, i) {return nfa.accept(q[i]);}).length;
        }
        if (',' === type)
        {
            return (q[1]+1 === input.length) && input[q[1]].accept(q[0]);
        }
    },

    reject: function(q) {
        var self = this, type = self.type, input = self.input;
        if (is_obj(type))
        {
            if (isset(type, 'min'))
            {
                return ((q[1] >= type.max) && input.accept(q[0])) || ((q[1] < type.min) && input.reject(q[0]));
            }
            else
            {
                return !q[0].length;
            }
        }
        if ('l' === type)
        {
            return 0 === q;
        }
        if ('^' === type)
        {
            return !q;
        }
        if ('$' === type)
        {
            return !q;
        }
        if ('!' === type)
        {
            return input.accept(q);
        }
        if ('|' === type)
        {
            return input.length === input.filter(function(nfa, i) {return nfa.reject(q[i]);}).length;
        }
        if (',' === type)
        {
            return input[q[1]].reject(q[0]);
        }
    },

    match: function(string, offset, q) {
        var self = this;
        var i = offset || 0,
            j = i,
            n = string.length;
        if (null == q) q = self.q0();
        for (;;)
        {
            if (j >= n)
            {
                ++i;
                if (i > n) break;
                j = i;
                q = self.q0();
            }
            if (0 === j) q = self.d(q, 0);
            q = self.d(q, string.charAt(j));
            if (n-1 === j) q = self.d(q, 1);
            if (self.accept(q))
            {
                return i; // matched
            }
            else if (self.reject(q))
            {
                j = n; // failed, try next
            }
            else
            {
                ++j; // continue
            }
        }
        return -1;
    }
};
Matchy.NFA = NFA;

// utils
var stdMath = Math,
    min = stdMath.min,
    max = stdMath.max,
    INF = Infinity,
    HAS = Object.prototype.hasOwnProperty,
    toString = Object.prototype.toString;

function isset(o, x)
{
    return HAS.call(o, x);
}
function is_number(x)
{
    return 'number' === typeof x;
}
function is_string(x)
{
    return '[object String]' === toString.call(x);
}
function is_array(x)
{
    return '[object Array]' === toString.call(x);
}
function is_obj(x)
{
    return '[object Object]' === toString.call(x);
}
function array_fill(i0, n, v)
{
    var array = new Array(n);
    for (var j=0,i=i0; j<n; ++i,++j) array[i] = v;
    return array;
}
function range(start, end, step)
{
    if (null == step) step = 1;
    var range = [];
    while (start <= end)
    {
        range.push(start);
        start += step;
    }
    return range;
}

// export it
return Matchy;
});

