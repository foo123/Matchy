<?php
/**
*  Matchy
*  String searching algorithms for PHP, JavaScript, Python
*
*  @version: 1.0.0
*  https://github.com/foo123/Matchy
*
**/
if (!class_exists('Matchy', false))
{
class Matchy
{
    const VERSION = "1.0.0";

    public function __construct()
    {
    }

    public function fsa($pattern, $string = null, $offset = 0)
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
        $suffix = array();
        $is_suffix = function($k, $q, $c) use ($pattern, $in_pattern, &$suffix) {
            if (!isset($in_pattern[$c])) return false;
            if (!isset($suffix[$k])) $suffix[$k] = array();
            if (!isset($suffix[$k][$q])) $suffix[$k][$q] = array();
            if (!isset($suffix[$k][$q][$c]))
            {
                $s1 = mb_substr($pattern, 0, $k, 'UTF-8');
                $s2 = mb_substr($pattern, 0, $q, 'UTF-8') . $c;
                $suffix[$k][$q][$c] = ($s1 === mb_substr($s2, -mb_strlen($s1, 'UTF-8'), null, 'UTF-8'));
            }
            return $suffix[$k][$q][$c];
        };

        $d = array_map(function() {return array();}, array_fill(0, $m+1, 0));
        $delta = function($q, $c) use (&$d, $pattern, $m, &$is_suffix) {
            if (isset($d[$q][$c])) return $d[$q][$c];
            $k = min($m, $q+1);
            while ((0 < $k) && !$is_suffix($k, $q, $c)) --$k;
            $d[$q][$c] = $k;
            return $k;
        };

        // matcher
        $matcher = function($string, $offset = 0) use ($m, &$delta) {
            $n = mb_strlen($string, 'UTF-8');
            if (0 > $offset) $offset += $n;
            if ((0 < $n) && (0 < $m) && ($n >= $offset+$m))
            {
                for ($q=0,$i=$offset; $i<$n; ++$i)
                {
                    $c = mb_substr($string, $i, 1, 'UTF-8');
                    $q = $delta($q, $c);
                    if ($m === $q) return $i-$m+1; // matched
                }
            }
            return -1;
        };
        return null === $string ? $matcher : $matcher($string, $offset);
    }

    public function rabinkarp($pattern, $string = null, $offset = 0)
    {
        // https://en.wikipedia.org/wiki/Rabin%E2%80%93Karp_algorithm
        // http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=BA184C94C16CB23D5FA7329E257E3713?doi=10.1.1.86.9502&rep=rep1&type=pdf

        $m = mb_strlen($pattern, 'UTF-8');
        $Q = 3989; // prime number so that 10*Q can fit in one WORD (i.e 2^16 bits)

        // matcher
        $matcher = function($string, $offset = 0) use ($pattern, $m, $Q) {
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
                        if (mb_substr($string, $i, $m, 'UTF-8') === $pattern) return $i; // matched
                    }
                    // update text hash for next char using Horner algorithm
                    if ($i < $n)
                    {
                        $sq = ($D*($sq - $h*$alphabet[mb_substr($string, $i, 1, 'UTF-8')]) + $alphabet[mb_substr($string, $i+$m, 1, 'UTF-8')]) % $Q;
                        if (0 > $sq) $sq += $Q;
                    }
                }
            }
            return -1;
        };
        return null === $string ? $matcher : $matcher($string, $offset);
    }

    public function knuthmorrispratt($pattern, $string = null, $offset = 0)
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
        $matcher = function($string, $offset = 0) use ($pattern, $m, $T) {
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
                        if ($k === $m) return $i - $k; // matched
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
            return -1;
        };
        return null === $string ? $matcher : $matcher($string, $offset);
    }

    public function boyermoore($pattern, $string = null, $offset = 0)
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
        $matcher = function($string, $offset = 0) use ($pattern, $m, $shift) {
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
                        return $i; // matched
                    }
                    else
                    {
                        $i += $shift[$j+1];
                    }
                }
            }
            return -1;
        };
        return null === $string ? $matcher : $matcher($string, $offset);
    }

    public function twoway($pattern, $string = null, $offset = 0)
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
        $matcher = function($string, $offset = 0) use ($pattern, $m, $small_period, $l, $p) {
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
                            if ($j <= $s) return $pos; // matched
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
                            if (0 === $j) return $pos; // matched
                            $pos += $q;
                        }
                    }
                }
            }
            return -1;
        };
        return null === $string ? $matcher : $matcher($string, $offset);
    }

    public function commentzwalter($pattern, $string = null, $offset = 0)
    {
        // https://en.wikipedia.org/wiki/Commentz-Walter_algorithm
        $non_matcher = function($string, $offset = 0) {
            return -1;
        };
        return null === $string ? $non_matcher : $non_matcher($string, $offset); // TODO
    }

    public function baezayatesgonnet($pattern, $string = null, $offset = 0)
    {
        // https://en.wikipedia.org/wiki/Bitap_algorithm
        $non_matcher = function($string, $offset = 0) {
            return -1;
        };
        return null === $string ? $non_matcher : $non_matcher($string, $offset); // TODO
    }

    public function ahocorasick($pattern, $string = null, $offset = 0)
    {
        // https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_algorithm
        $non_matcher = function($string, $offset = 0) {
            return -1;
        };
        return null === $string ? $non_matcher : $non_matcher($string, $offset); // TODO
    }
}
}
