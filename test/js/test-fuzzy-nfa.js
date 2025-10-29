"use strict";

const Matchy = require('../../src/js/Matchy.js');
const echo = console.log;

function NFA(input, type = 'l')
{
    return new Matchy.NFA(input, type);
}
function create_string(alphabet, n)
{
    let s = '';
    for (let i=0; i<n; ++i)
    {
        s += alphabet[Math.round(Math.random()*(alphabet.length-1))];
    }
    return s;
}
function test_case(nfa, pattern, string, match, errors)
{
    const found = nfa.match(string);
    echo('fuzzynfa("'+pattern+'", "'+string+'") = '+found[0]+', errors '+found[1]+' (expected '+match+', errors '+errors+')');
}
function test()
{
    let test;
    test = {pattern:'^(a+)(b+)$', nfa:NFA(NFA([NFA('', '^'), NFA(NFA('a'), '+'), NFA(NFA('b'), '+'), NFA('', '$')], ','), {total_errors:2,word_level:false})};
    test_case(test.nfa, test.pattern, "aaabbbbb", 0, 0); // 0 errors
    test_case(test.nfa, test.pattern, "ababbbbb", 0, 1); // 1 errors
    test_case(test.nfa, test.pattern, "abababbb", 0, 2); // 2 errors
    test_case(test.nfa, test.pattern, "abababab", -1, 3); // 3 errors
    test_case(test.nfa, test.pattern, "aabababbbb", 0, 2); // 2 errors
    test_case(test.nfa, test.pattern, "baabaaabbb", 0, 2); // 2 errors
    
    test = {pattern:'(aaa)(bbb)', nfa:NFA(NFA([NFA('aaa'), NFA('bbb')], ','), {total_errors:1,word_level:true})};
    test_case(test.nfa, test.pattern, "aaabbb", 0, 0); // 0 errors
    test_case(test.nfa, test.pattern, "bbb", 0, 1); // 1 errors, deletion
    test_case(test.nfa, test.pattern, "cbbb", 1, 1); // 1 errors, deletion or substitution
    test_case(test.nfa, test.pattern, "aacbbb", 3, 1); // 1 errors, deletion or substitution
    test_case(test.nfa, test.pattern, "aaacbbb", 0, 1); // 1 errors, insertion
   
    test = {pattern:'^(aaa)(bbb)$', nfa:NFA(NFA([NFA('', '^'), NFA('aaa'), NFA('bbb'), NFA('', '$')], ','), {total_errors:1,word_level:true,transpositions:true})};
    test_case(test.nfa, test.pattern, "aaabbb", 0, 0); // 0 errors
    test_case(test.nfa, test.pattern, "bbbaaa", 0, 1); // 1 errors, transposition
    test_case(test.nfa, test.pattern, "bbb", 0, 1); // 1 errors, deletion
    test_case(test.nfa, test.pattern, "cbbb", 0, 1); // 1 errors, substitution
    test_case(test.nfa, test.pattern, "aacbbb", 0, 1); // 1 errors, substitution
    test_case(test.nfa, test.pattern, "aaacbbb", 0, 1); // 1 errors, insertion
}

test();
