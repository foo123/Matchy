"use strict";

const Matchy = require('../../src/js/Matchy.js');
const echo = console.log;
const NFA = Matchy.NFA;

function create_string(alphabet, n)
{
    let s = '';
    for (let i=0; i<n; ++i)
    {
        s += alphabet[Math.round(Math.random()*(alphabet.length-1))];
    }
    return s;
}
function test_case(nfa, pattern, string, match)
{
    const found = nfa.match(string);
    echo('fuzzynfa("'+pattern+'", "'+string+'") = '+found+' ('+match+')');
}
function test()
{
    // some special test cases
    let test;
    test = {pattern:'^(a+)(b+)$', nfa:NFA(NFA([NFA('', '^'), NFA(NFA('a'), '+'), NFA(NFA('b'), '+'), NFA('', '$')], ','), {total_errors:2})};
    test_case(test.nfa, test.pattern, "aaabbbbb", 0); // 0 errors
    test_case(test.nfa, test.pattern, "ababbbbb", 0); // 1 errors
    test_case(test.nfa, test.pattern, "abababbb", 0); // 2 errors
    test_case(test.nfa, test.pattern, "abababab", -1); // 3 errors
    test_case(test.nfa, test.pattern, "aabababbbb", 0); // 2 errors
    test_case(test.nfa, test.pattern, "baabaaabbb", 0); // 2 errors
}

test();
