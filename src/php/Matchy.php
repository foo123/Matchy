<?php
/**
*  Matchy
*  Exact and fuzzy string searching algorithms for PHP, JavaScript, Python
*
*  @version: 4.0.0
*  https://github.com/foo123/Matchy
*
**/
if (!class_exists('Matchy', false))
{
class Matchy
{
    const VERSION = "4.0.0";

    public function __construct()
    {
    }

    public function fsa($pattern, $string = null, $offset = 0, $return_match = false)
    {
        // https://en.wikipedia.org/wiki/Finite-state_machine
        // https://en.wikipedia.org/wiki/String-searching_algorithm
        // https://euroinformatica.ro/documentation/programming/!!!Algorithms_CORMEN!!!/DDU0214.html

        // construct transition matrix
        $m = mb_strlen($pattern, 'UTF-8');
        $in_pattern = array();
        for ($i=0; $i<$m; ++$i)
        {
            $in_pattern[mb_substr($pattern, $i, 1, 'UTF-8')] = 1;
        }
        $is_suffix = function($k, $q, $c) use ($pattern) {
            $s1 = mb_substr($pattern, 0, $k, 'UTF-8');
            $s2 = mb_substr($pattern, 0, $q, 'UTF-8') . $c;
            return $s1 === mb_substr($s2, -mb_strlen($s1, 'UTF-8'), null, 'UTF-8');
        };

        $d = array_map(function() {return array();}, array_fill(0, $m+1, 0));
        $delta = function($q, $c) use (&$d, $pattern, $m, $in_pattern, &$is_suffix) {
            if (!isset($in_pattern[$c])) return 0;
            if (!isset($d[$q][$c]))
            {
                $k = min($m, $q+1);
                while ((0 < $k) && !$is_suffix($k, $q, $c)) --$k;
                $d[$q][$c] = $k;
            }
            return $d[$q][$c];
        };

        // matcher
        $matcher = function($string, $offset = 0, $return_match = false) use ($pattern, $m, &$delta) {
            $n = mb_strlen($string, 'UTF-8');
            if (0 > $offset) $offset += $n;
            if ((0 < $n) && (0 < $m) && ($n >= $offset+$m))
            {
                for ($q=0,$i=$offset; $i<$n; ++$i)
                {
                    $c = mb_substr($string, $i, 1, 'UTF-8');
                    $q = $delta($q, $c);
                    if ($m === $q) return $return_match ? $pattern : ($i-$m+1); // matched
                }
            }
            return $return_match ? null : -1;
        };
        return null === $string ? $matcher : $matcher($string, $offset, $return_match);
    }

    public function rabinkarp($pattern, $string = null, $offset = 0, $return_match = false)
    {
        // https://en.wikipedia.org/wiki/Rabin%E2%80%93Karp_algorithm
        // http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=BA184C94C16CB23D5FA7329E257E3713?doi=10.1.1.86.9502&rep=rep1&type=pdf

        $m = mb_strlen($pattern, 'UTF-8');
        $Q = 3989; // prime number so that 10*Q can fit in one WORD (i.e 2^16 bits)

        // matcher
        $matcher = function($string, $offset = 0, $return_match = false) use ($pattern, $m, $Q) {
            $n = mb_strlen($string, 'UTF-8');
            if (0 > $offset) $offset += $n;
            if ((0 < $n) && (0 < $m) && ($n >= $offset+$m))
            {
                $alphabet = array();
                $D = 0;
                for ($i=0; $i<$m; ++$i)
                {
                    $c = mb_substr($pattern, $i, 1, 'UTF-8');
                    if (!isset($alphabet[$c])) $alphabet[$c] = $D++;
                }
                for ($i=$offset; $i<$n; ++$i)
                {
                    $c = mb_substr($string, $i, 1, 'UTF-8');
                    if (!isset($alphabet[$c])) $alphabet[$c] = $D++;
                }
                $h = 1;
                $pq = $alphabet[mb_substr($pattern, 0, 1, 'UTF-8')] % $Q;
                $sq = $alphabet[mb_substr($string, $offset, 1, 'UTF-8')] % $Q;
                for ($i=1; $i<$m; ++$i)
                {
                    $h  = ($h*$D) % $Q;
                    $pq = ($pq*$D + $alphabet[mb_substr($pattern, $i, 1, 'UTF-8')]) % $Q;
                    $sq = ($sq*$D + $alphabet[mb_substr($string, $offset+$i, 1, 'UTF-8')]) % $Q;
                }
                $n = $n-$m;
                for ($i=$offset; $i<=$n; ++$i)
                {
                    // pq, sq, D-base arithmetic code, used as a quick "hash" test
                    // worst case: many hash collisions -> "naive" matching
                    if ($pq === $sq)
                    {
                        if (mb_substr($string, $i, $m, 'UTF-8') === $pattern) return $return_match ? $pattern : $i; // matched
                    }
                    // update text hash for next char using Horner algorithm
                    if ($i < $n)
                    {
                        $sq = ($D*($sq - $h*$alphabet[mb_substr($string, $i, 1, 'UTF-8')]) + $alphabet[mb_substr($string, $i+$m, 1, 'UTF-8')]) % $Q;
                        if (0 > $sq) $sq += $Q;
                    }
                }
            }
            return $return_match ? null : -1;
        };
        return null === $string ? $matcher : $matcher($string, $offset, $return_match);
    }

    public function knuthmorrispratt($pattern, $string = null, $offset = 0, $return_match = false)
    {
        // https://en.wikipedia.org/wiki/Knuth%E2%80%93Morris%E2%80%93Pratt_algorithm
        // http://www.eecs.ucf.edu/~shzhang/Combio09/kmp.pdf

        // pre-processing
        $m = mb_strlen($pattern, 'UTF-8');
        $T = array_fill(0, $m, 0);
        $T[0] = -1;
        $k = 1;
        $i = 0;
        while ($k < $m)
        {
            $c = mb_substr($pattern, $k, 1, 'UTF-8');
            if ($c === mb_substr($pattern, $i, 1, 'UTF-8'))
            {
                $T[$k] = $T[$i];
            }
            else
            {
                $T[$k] = $i;
                while ((0 <= $i) && ($c !== mb_substr($pattern, $i, 1, 'UTF-8'))) $i = $T[$i];
            }
            ++$k;
            ++$i;
        }

        // matcher
        $matcher = function($string, $offset = 0, $return_match = false) use ($pattern, $m, $T) {
            $n = mb_strlen($string, 'UTF-8');
            if (0 > $offset) $offset += $n;
            if ((0 < $n) && (0 < $m) && ($n >= $offset+$m))
            {
                $k = 0;
                $i = $offset;
                while ($i < $n)
                {
                    if (mb_substr($pattern, $k, 1, 'UTF-8') === mb_substr($string, $i, 1, 'UTF-8'))
                    {
                        ++$i;
                        ++$k;
                        if ($k === $m) return $return_match ? $pattern : ($i - $k); // matched
                    }
                    else
                    {
                        $k = $T[$k];
                        if ($k < 0)
                        {
                            ++$i;
                            ++$k;
                        }
                    }
                }
            }
            return $return_match ? null : -1;
        };
        return null === $string ? $matcher : $matcher($string, $offset, $return_match);
    }

    public function boyermoore($pattern, $string = null, $offset = 0, $return_match = false)
    {
        // https://en.wikipedia.org/wiki/Boyer%E2%80%93Moore_string-search_algorithm
        // http://www.cs.utexas.edu/~moore/publications/fstrpos.pdf

        // pre-processing
        $m = mb_strlen($pattern, 'UTF-8');
        $suffix = array_fill(0, $m+1, 0);
        $shift = array_fill(0, $m+1, 0);
        $i = $m;
        $j = $m+1;
        $suffix[$i] = $j;
        while (0 < $i)
        {
            while (($j <= $m) && (mb_substr($pattern, $i-1, 1, 'UTF-8') !== mb_substr($pattern, $j-1, 1, 'UTF-8')))
            {
                if (0 === $shift[$j]) $shift[$j] = $j-$i;
                $j = $suffix[$j];
            }
            --$i;
            --$j;
            $suffix[$i] = $j;
        }
        $j = $suffix[0];
        for ($i=0; $i<$m; ++$i)
        {
            if (0 === $shift[$i])
            {
                $shift[$i] = $j;
                if ($i === $j) $j = $suffix[$j];
            }
        }

        // matcher
        $matcher = function($string, $offset = 0, $return_match = false) use ($pattern, $m, $shift) {
            $n = mb_strlen($string, 'UTF-8');
            if (0 > $offset) $offset += $n;
            if ((0 < $n) && (0 < $m) && ($n >= $offset+$m))
            {
                $i = $offset;
                while ($i + $m <= $n)
                {
                    $j = $m - 1;
                    while ((0 <= $j) && (mb_substr($pattern, $j, 1, 'UTF-8') === mb_substr($string, $i+$j, 1, 'UTF-8')))
                    {
                        --$j;
                    }
                    if (0 > $j)
                    {
                        return $return_match ? $pattern : $i; // matched
                    }
                    else
                    {
                        $i += $shift[$j+1];
                    }
                }
            }
            return $return_match ? null : -1;
        };
        return null === $string ? $matcher : $matcher($string, $offset, $return_match);
    }

    public function twoway($pattern, $string = null, $offset = 0, $return_match = false)
    {
        // https://en.wikipedia.org/wiki/Two-way_string-matching_algorithm
        // http://monge.univ-mlv.fr/~mac/Articles-PDF/CP-1991-jacm.pdf

        // pre-processing
        $m = mb_strlen($pattern, 'UTF-8');
        $maximal_suffix = function($s, $n, $dir) {
            // assumes 1-indexed strings instead of 0-indexed
            $p = 1;  //  currently known period.
            $k = 1;  // index for period testing, 0 < k <= p.
            $j = 1;  // index for maxsuf testing. greater than maxs.
            $i = 0;  // the proposed starting index of maxsuf
            while ($j + $k <= $n)
            {
                $a = mb_substr($s, $j + $k - 1, 1, 'UTF-8');
                $b = mb_substr($s, $i + $k - 1, 1, 'UTF-8');
                $cmp = $dir*($a > $b ? 1 : ($a < $b ? -1 : 0));
                if (0 > $cmp)
                {
                    // Suffix is smaller. Period is the entire prefix so far.
                    $j += $k;
                    $k = 1;
                    $p = $j - $i;
                }
                elseif (0 < $cmp)
                {
                    // Suffix is larger. Start over from here.
                    $i = $j;
                    $j = $i + 1;
                    $k = 1;
                    $p = 1;
                }
                else
                {
                    // They are the same - we should go on.
                    if ($k === $p)
                    {
                        // We are done checking this stretch of p. reset k.
                        $j += $p;
                        $k = 1;
                    }
                    else
                    {
                        ++$k;
                    }
                }
            }
            return array($i, $p);
        };
        // critical factorization
        $s1 = $maximal_suffix($pattern, $m,  1);
        $s2 = $maximal_suffix($pattern, $m, -1);
        if ($s1[0] >= $s2[0])
        {
            list($l, $p) = $s1;
        }
        else
        {
            list($l, $p) = $s2;
        }
        // small period
        $small_period = (2*$l < $m) && (mb_substr($pattern, 0, $l, 'UTF-8') === mb_substr(mb_substr($pattern, $l, $p, 'UTF-8'), -$l, null, 'UTF-8'));

        // matcher
        $matcher = function($string, $offset = 0, $return_match = false) use ($pattern, $m, $small_period, $l, $p) {
            $n = mb_strlen($string, 'UTF-8');
            if (0 > $offset) $offset += $n;
            if ((0 < $n) && (0 < $m) && ($n >= $offset+$m))
            {
                // assumes 1-indexed strings instead of 0-indexed
                $pos = $offset;
                if ($small_period)
                {
                    // small period
                    $s = 0;
                    while ($pos + $m <= $n)
                    {
                        $i = max($l, $s) + 1;
                        while (($i <= $m) && (mb_substr($pattern, $i-1, 1, 'UTF-8') === mb_substr($string, $pos+$i-1, 1, 'UTF-8')))
                        {
                            ++$i;
                        }
                        if ($i <= $m)
                        {
                            $pos += max($i-$l, $s-$p+1);
                            $s = 0;
                        }
                        else
                        {
                            $j = $l;
                            while (($j > $s) && (mb_substr($pattern, $j-1, 1, 'UTF-8') === mb_substr($string, $pos+$j-1, 1, 'UTF-8')))
                            {
                                --$j;
                            }
                            if ($j <= $s) return $return_match ? $pattern : $pos; // matched
                            $pos += $p;
                            $s = $m-$p;
                        }
                    }
                }
                else
                {
                    // large period
                    $q = max($l, $m-$l) + 1;
                    while ($pos + $m <= $n)
                    {
                        $i = $l + 1;
                        while (($i <= $m) && (mb_substr($pattern, $i-1, 1, 'UTF-8') === mb_substr($string, $pos+$i-1, 1, 'UTF-8')))
                        {
                            ++$i;
                        }
                        if ($i <= $m)
                        {
                            $pos += $i-$l;
                        }
                        else
                        {
                            $j = $l;
                            while (($j > 0) && (mb_substr($pattern, $j-1, 1, 'UTF-8') === mb_substr($string, $pos+$j-1, 1, 'UTF-8')))
                            {
                                --$j;
                            }
                            if (0 === $j) return $return_match ? $pattern : $pos; // matched
                            $pos += $q;
                        }
                    }
                }
            }
            return $return_match ? null : -1;
        };
        return null === $string ? $matcher : $matcher($string, $offset, $return_match);
    }

    public function commentzwalter($pattern, $string = null, $offset = 0, $return_match = false)
    {
        // https://en.wikipedia.org/wiki/Commentz-Walter_algorithm
        $non_matcher = function($string, $offset = 0, $return_match = false) {
            return $return_match ? null : -1;
        };
        return null === $string ? $non_matcher : $non_matcher($string, $offset, $return_match); // TODO
    }

    public function baezayatesgonnet($pattern, $string = null, $offset = 0, $return_match = false)
    {
        // https://en.wikipedia.org/wiki/Bitap_algorithm
        $non_matcher = function($string, $offset = 0, $return_match = false) {
            return $return_match ? null : -1;
        };
        return null === $string ? $non_matcher : $non_matcher($string, $offset, $return_match); // TODO
    }

    public function ahocorasick($pattern, $string = null, $offset = 0, $return_match = false)
    {
        // https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_algorithm
        $non_matcher = function($string, $offset = 0, $return_match = false) {
            return $return_match ? null : -1;
        };
        return null === $string ? $non_matcher : $non_matcher($string, $offset, $return_match); // TODO
    }
}
// non-deterministic finite automaton
class MatchyNFA
{
    public static function makeFuzzy($nfa, $max_errors, $at_word_level = false, $transpositions = false)
    {
        if ($at_word_level)
        {
            if ($nfa instanceof MatchyNFA)
            {
                if (is_array($nfa->type) && (isset($nfa->type['errors']) || isset($nfa->type['total_errors'])))
                {
                    // leave as is
                }
                elseif (',' === $nfa->type)
                {
                    // approximate at word level
                    $nfa = new MatchyNFA(MatchyNFA::makeFuzzy($nfa->input, $max_errors, $at_word_level, $transpositions), $transpositions ? ',,,' : ',,');
                }
                elseif (is_array($nfa->input))
                {
                    // recurse
                    $nfa->input = array_map(function($input) use ($max_errors, $at_word_level, $transpositions) {
                        return MatchyNFA::makeFuzzy($input, $max_errors, $at_word_level, $transpositions);
                    }, $nfa->input);
                }
                elseif ($nfa->input instanceof MatchyNFA)
                {
                    // recurse
                    $nfa->input = MatchyNFA::makeFuzzy($nfa->input, $max_errors, $at_word_level, $transpositions);
                }
            }
        }
        else // at char level
        {
            if ($nfa instanceof MatchyNFA)
            {
                if (is_array($nfa->type) && (isset($nfa->type['errors']) || isset($nfa->type['total_errors'])))
                {
                    // leave as is
                }
                elseif ('l' === $nfa->type)
                {
                    // change literal match to approximate match
                    $nfa = new MatchyNFA($nfa->input, array('errors'=>$max_errors, 'transpositions'=>$transpositions));
                }
                elseif (is_array($nfa->input))
                {
                    // recurse
                    $nfa->input = array_map(function($input) use ($max_errors, $at_word_level, $transpositions) {
                        return MatchyNFA::makeFuzzy($input, $max_errors, $at_word_level, $transpositions);
                    }, $nfa->input);
                }
                elseif ($nfa->input instanceof MatchyNFA)
                {
                    // recurse
                    $nfa->input = MatchyNFA::makeFuzzy($nfa->input, $max_errors, $at_word_level, $transpositions);
                }
            }
            elseif (is_string($nfa))
            {
                // change literal string to approximate match
                $nfa = new MatchyNFA($nfa, array('errors'=>$max_errors, 'transpositions'=>$transpositions));
            }
        }
        return $nfa;
    }

    public $input = null;
    public $type = null;

    public function __construct($input, $type = 'l')
    {
        if ('?' === $type) $type = array('min'=>0, 'max'=>1);
        if ('*' === $type) $type = array('min'=>0, 'max'=>INF);
        if ('+' === $type) $type = array('min'=>1, 'max'=>INF);
        $this->type = $type;
        $this->input = is_array($type) && isset($type['total_errors']) ? MatchyNFA::makeFuzzy($input, $type['total_errors'], !empty($type['word_level']), !empty($type['transpositions'])) : $input;
    }

    public function q0()
    {
        $type = $this->type;
        $input = $this->input;
        if (is_array($type))
        {
            if (isset($type['min']))
            {
                $q = array($input->q0(), 0);
            }
            elseif (isset($type['errors']))
            {
                $k = min($type['errors'], mb_strlen($input, 'UTF-8'));
                $q = !empty($type['transpositions']) ? array(
                    range(0, $k, 1),
                    range(0, $k, 1),
                    array(),
                    array(),
                    ''
                ) : array(
                    range(0, $k, 1),
                    range(0, $k, 1)
                );
            }
            else
            {
                $q = $input->q0();
            }
        }
        if ('l' === $type)
        {
            $q = 0;
        }
        if ('^' === $type)
        {
            $q = false;
        }
        if ('$' === $type)
        {
            $q = false;
        }
        if ('!' === $type)
        {
            $q = $input->q0();
        }
        if ('|' === $type)
        {
            $q = array_map(function($nfa) {
                return $nfa->q0();
            }, $input);
        }
        if (',' === $type)
        {
            $q = array(array($input[0]->q0(), 0, 0));
        }
        if ((',,' === $type) || (',,,' === $type))
        {
            $q = array();
            for ($i=0,$n=count($input); $i<$n; ++$i)
            {
                // push for insertion, deletion, substitution, transposition
                $q[] = array($input[$i]->q0(), 0, $i, $i);
            }
        }
        return array('q'=>$q, 'e'=>0); // keep track of errors
    }

    public function d($qe, $c)
    {
        $type = $this->type;
        $input = $this->input;
        $q = $qe['q'];
        $e = $qe['e'];
        if (is_array($type))
        {
            if (isset($type['min']))
            {
                $e0 = $q[0]['e'];
                $qe = $input->d($q[0], $c);
                $q = array($qe, $q[1]);
                $e += $q[0]['e'] - $e0;
                if ($input->accept($qe)) $q = array($input->q0(), $q[1]+1);
            }
            elseif (isset($type['errors']))
            {
                if (is_int($c))
                {
                    //$q = $q;
                }
                else
                {
                    $transpositions = !empty($type['transpositions']);
                    $w = $input;
                    $n = mb_strlen($w, 'UTF-8');
                    $k = min($type['errors'], $n);
                    $min_e = $k+1;
                    $index = $q[0];
                    $value = $q[1];
                    $new_index = array();
                    $new_value = array();
                    $m = count($index);
                    $prev_i = -1;
                    $prev_v = 0;
                    $next_i = -1;
                    if ($transpositions)
                    {
                        $index_2 = $q[2];
                        $value_2 = $q[3];
                        $cp = $q[4];
                        $m2 = count($index_2);
                        $j2 = 0;
                    }
                    if ((0 < $m) && (0 === $index[0]) && ($value[0] < $k))
                    {
                        $i = 0;
                        $v = $value[0] + 1;
                        $prev_i = $i;
                        $prev_v = $v;
                        $new_index[] = $i;
                        $new_value[] = $v;
                        $min_e = min($min_e, $v);
                    }
                    foreach ($index as $j => $i)
                    {
                        if ($i >= $n) break;
                        $d = mb_substr($w, $i, 1, 'UTF-8') === $c ? 0 : 1;
                        $v = $value[$j] + $d; // L[i,ii] = L[i-1,ii-1] + d
                        $next_i = $j+1 < $m ? $index[$j+1] : -1;
                        ++$i;
                        if ($i-1 === $prev_i)
                        {
                            $v = min($v, $prev_v + 1); // L[i,ii] = min(L[i,ii], L[i-1,ii] + 1)
                        }
                        if ($i === $next_i)
                        {
                            $v = min($v, $value[$j+1] + 1); // L[i,ii] = min(L[i,ii], L[i,ii-1] + 1)
                        }
                        if ($transpositions && ($cp === mb_substr($w, $i-1, 1, 'UTF-8')) && ($c === mb_substr($w, $i-2, 1, 'UTF-8')))
                        {
                            while (($j2 < $m2) && ($index_2[$j2] < $i-2)) ++$j2;
                            if (($j2 < $m2) && ($i-2 === $index_2[$j2]))
                            {
                                $v = min($v, $value_2[$j2] + $d); // L[i,ii] = min(L[i,ii], L[i-2,ii-2] + d)
                                ++$j2;
                            }
                        }
                        if ($v <= $k)
                        {
                            $min_e = min($min_e, $v);
                            $prev_i = $i;
                            $prev_v = $v;
                            $new_index[] = $i;
                            $new_value[] = $v;
                        }
                    }
                    $q = $transpositions ? array(
                        $new_index,
                        $new_value,
                        $index,
                        $value,
                        $c
                    ) : array(
                        $new_index,
                        $new_value
                    );
                    $e = $min_e;
                }
            }
            else
            {
                $q = $input->d($q, $c);
                $e = $q['e'];
            }
        }
        if ('l' === $type)
        {
            if (is_int($c))
            {
                //$q = $q;
            }
            else
            {
                $q = $q < mb_strlen($input, 'UTF-8') ? ($c === mb_substr($input, $q, 1, 'UTF-8') ? $q+1 : 0) : 0;
            }
        }
        if ('^' === $type)
        {
            $q = is_int($c) && (0 === $c);
            //$e = $q ? 0 : 1;
        }
        if ('$' === $type)
        {
            $q = is_int($c) && (1 === $c);
            //$e = $q ? 0 : 1;
        }
        if ('!' === $type)
        {
            $q = $input->d($q, $c);
            $e = $q['e'];
        }
        if ('|' === $type)
        {
            $e0 = $e;
            $e = INF;
            $q = array_map(function($i) use ($input, $q, $c) {
                return $input[$i]->d($q[$i], $c);
            }, array_keys($input));
            foreach ($q as $i => $qi)
            {
                if (!$input[$i]->reject($qi)) $e = min($e, $qi['e']);
            }
            if (!is_finite($e)) $e = $e0+1;
        }
        if (',' === $type)
        {
            $n = count($input);
            $qq = $q;
            $q = array();
            foreach ($qq as $qi)
            {
                $i = $qi[1];
                $ei = $qi[2];
                $e0 = $qi[0]['e'];
                if ($input[$i]->accept($qi[0]))
                {
                    if ($i+1 < $n)
                    {
                        $q1 = $input[$i]->d($qi[0], $c);
                        if (!$input[$i]->reject($q1))
                        {
                            $qi = array($q1, $i, $ei);
                            $q[] = $qi;
                        }
                        $q1 = $input[$i+1]->d($input[$i+1]->q0(), $c);
                        if (!$input[$i+1]->reject($q1))
                        {
                            $qi = array($q1, $i+1, $ei+$e0);
                            $q[] = $qi;
                        }
                    }
                    else
                    {
                        $q[] = $qi;
                    }
                }
                else
                {
                    $qi = array($input[$i]->d($qi[0], $c), $i, $ei);
                    $q[] = $qi;
                }
            }
            $last_i = 0;
            foreach ($q as $qi)
            {
                $i = $qi[1];
                if ($i > $last_i)
                {
                    $last_i = $i;
                    $e = $qi[0]['e'] + $qi[2];
                }
                elseif ($i === $last_i)
                {
                    $e = min($e, $qi[0]['e'] + $qi[2]);
                }
            }
        }
        if ((',,' === $type) || (',,,' === $type))
        {
            $n = count($input);
            $qq = $q;
            $q = array();
            foreach ($qq as $qi)
            {
                $i = $qi[1];
                $j = $qi[2];
                $next_i = $i+1 < $n ? $i+1 : $i;
                $ei = $qi[3];
                if ($input[$j]->accept($qi[0]))
                {
                    $q1 = $input[$j]->d($qi[0], $c);
                    if (!$input[$j]->reject($q1))
                    {
                        $qi = array($q1, $i, $j, $ei);
                        $q[] = $qi;
                    }
                    else
                    {
                        $q[] = $qi;
                    }
                    // push next for insertion, deletion, substitution
                    for ($k=1; $j+$k<$n; ++$k)
                    {
                        $q1 = $input[$j+$k]->d($input[$j+$k]->q0(), $c);
                        if (!$input[$j+$k]->reject($q1))
                        {
                            $qi = array($q1, $next_i, $j+$k, $ei+abs($j+$k-$next_i));
                            $q[] = $qi;
                        }
                    }
                    // push prev for transposition
                    if ((',,,' === $type) && (0 < $j) && ($i+1 === $j))
                    {
                        $q1 = $input[$j-1]->d($input[$j-1]->q0(), $c);
                        if (!$input[$j-1]->reject($q1))
                        {
                            $qi = array($q1, $next_i, $j-1, $ei-1+1); // carries the previous error of deletion
                            $q[] = $qi;
                        }
                    }
                }
                else
                {
                    $q1 = $input[$j]->d($qi[0], $c);
                    if ($input[$j]->reject($q1))
                    {
                        // push next for insertion, substitution, deletion
                        for ($k=0; $j+$k<$n; ++$k)
                        {
                            $q1 = $input[$j+$k]->d($input[$j+$k]->q0(), $c);
                            if (!$input[$j+$k]->reject($q1))
                            {
                                $qi = array($q1, $next_i, $j+$k, $ei+abs($j+$k-$next_i));
                                $q[] = $qi;
                            }
                        }
                    }
                    else
                    {
                        $qi = array($q1, $i, $j, $ei);
                        $q[] = $qi;
                    }
                }
            }
            $e0 = $e;
            $e = INF;
            foreach ($q as $qi)
            {
                if (($qi[1]+1 === $n) && $input[$qi[2]]->accept($qi[0]))
                {
                    $e = min($e, $qi[3]);
                }
            }
            if (!is_finite($e)) $e = $e0+1;
        }
        return array('q'=>$q, 'e'=>$e); // keep track of errors
    }

    public function accept($qe)
    {
        $type = $this->type;
        $input = $this->input;
        $q = $qe['q'];
        $e = $qe['e'];
        if (is_array($type))
        {
            if (isset($type['min']))
            {
                return ($type['min'] <= $q[1]) && ($q[1] <= $type['max']);
            }
            elseif (isset($type['errors']))
            {
                return (0 < count($q[0])) && ($q[0][count($q[0])-1] === mb_strlen($input, 'UTF-8'));
            }
            else
            {
                return $input->accept_with_errors($q, $type['total_errors']);
            }
        }
        if ('l' === $type)
        {
            return $q === mb_strlen($input, 'UTF-8');
        }
        if ('^' === $type)
        {
            return $q;
        }
        if ('$' === $type)
        {
            return $q;
        }
        if ('!' === $type)
        {
            return $input->reject($q);
        }
        if ('|' === $type)
        {
            return 0 < count(array_filter(array_keys($input), function($i) use ($input, $q) {
                return $input[$i]->accept($q[$i]);
            }));
        }
        if (',' === $type)
        {
            $n = count($input);
            return 0 < count(array_filter($q, function($qi) use ($n, $input) {
                return ($qi[1]+1 === $n) && $input[$qi[1]]->accept($qi[0]);
            }));
        }
        if ((',,' === $type) || (',,,' === $type))
        {
            $n = count($input);
            return 0 < count(array_filter($q, function($qi) use ($n, $input) {
                return ($qi[1]+1 === $n) && $input[$qi[2]]->accept($qi[0]);
            }));
        }
    }

    public function reject($qe)
    {
        $type = $this->type;
        $input = $this->input;
        $q = $qe['q'];
        $e = $qe['e'];
        if (is_array($type))
        {
            if (isset($type['min']))
            {
                return (($q[1] >= $type['max']) && $input->accept($q[0])) || (($q[1] < $type['min']) && $input->reject($q[0]));
            }
            elseif (isset($type['errors']))
            {
                return empty($q[0]);
            }
            else
            {
                return $input->reject_with_errors($q, $type['total_errors']);
            }
        }
        if ('l' === $type)
        {
            return 0 === $q;
        }
        if ('^' === $type)
        {
            return !$q;
        }
        if ('$' === $type)
        {
            return !$q;
        }
        if ('!' === $type)
        {
            return $input->accept($q);
        }
        if ('|' === $type)
        {
            return count($input) === count(array_filter(array_keys($input), function($i) use ($input, $q) {
                return $input[$i]->reject($q[$i]);
            }));
        }
        if (',' === $type)
        {
            return count($q) === count(array_filter($q, function($qi) use ($input) {
                return $input[$qi[1]]->reject($qi[0]);
            }));
        }
        if ((',,' === $type) || (',,,' === $type))
        {
            return count($q) === count(array_filter($q, function($qi) use ($input) {
                return $input[$qi[2]]->reject($qi[0]);
            }));
        }
    }

    public function accept_with_errors($qe, $max_errors = null)
    {
        if (!$this->accept($qe)) return false;
        if (null === $max_errors) return true;
        $type = $this->type;
        $input = $this->input;
        $q = $qe['q'];
        $e = $qe['e'];
        if (is_array($type) || (',' !== substr($type, 0, 1)))
        {
            return $e <= $max_errors;
        }
        if (',' === $type)
        {
            $n = count($input);
            return 0 < count(array_filter($q, function($qi) use ($n, $input, $max_errors) {
                return ($qi[1]+1 === $n) && ($qi[0]['e'] + $qi[2] <= $max_errors) && $input[$qi[1]]->accept($qi[0]);
            }));
        }
        if ((',,' === $type) || (',,,' === $type))
        {
            $n = count($input);
            return 0 < count(array_filter($q, function($qi) use ($n, $input, $max_errors) {
                return ($qi[1]+1 === $n) && ($qi[3] <= $max_errors) && $input[$qi[2]]->accept($qi[0]);
            }));
        }
    }

    public function reject_with_errors($qe, $max_errors = null)
    {
        if ($this->reject($qe)) return true;
        if (null === $max_errors) return false;
        $type = $this->type;
        $input = $this->input;
        $q = $qe['q'];
        $e = $qe['e'];
        if (is_array($type) || (',' !== substr($type, 0, 1)))
        {
            return $e > $max_errors;
        }
        if (',' === $type)
        {
            return count($q) === count(array_filter($q, function($qi) use ($input, $max_errors) {
                return ($qi[0]['e'] + $qi[2] > $max_errors) || $input[$qi[1]]->reject($qi[0]);
            }));
        }
        if ((',,' === $type) || (',,,' === $type))
        {
            return count($q) === count(array_filter($q, function($qi) use ($input, $max_errors) {
                return ($qi[3] > $max_errors) || $input[$qi[2]]->reject($qi[0]);
            }));
        }
    }

    public function get_errors($qe)
    {
        return $qe['e'];
    }

    public function match($string, $offset = 0, $return_match = false, $q = null)
    {
        $i = $offset;
        $j = $i;
        $n = mb_strlen($string, 'UTF-8');
        $c = '';
        $e = 0;
        if (null === $q) $q = $this->q0();
        for (;;)
        {
            if ($j >= $n)
            {
                ++$i;
                if ($i >= $n) break;
                $j = $i;
                $q = $this->q0();
                $c = mb_substr($string, $j-1, 1, 'UTF-8');
            }
            $prevc = $c;
            $c = mb_substr($string, $j, 1, 'UTF-8');
            if ((0 === $j) || ("\n" === $prevc)) $q = $this->d($q, 0);
            $q = $this->d($q, $c);
            if (($n-1 === $j) || ("\n" === $c)) $q = $this->d($q, 1);
            if ($this->accept($q))
            {
                $e = $this->get_errors($q);
                return array($return_match ? mb_substr($string, $i, $j-$i+1, 'UTF-8') : $i, $e); // matched
            }
            elseif ($this->reject($q))
            {
                $e = max($e, $this->get_errors($q));
                $j = $n; // failed, try next
            }
            else
            {
                ++$j; // continue
            }
        }
        return array($return_match ? null : -1, $e);
    }
}
}
