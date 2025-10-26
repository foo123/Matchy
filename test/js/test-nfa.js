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
function test_case(nfa, pattern, string, offset)
{
    const found = nfa.match(string, offset || 0);
    echo('nfa("'+pattern+'", "'+string+'", '+(offset||0)+') = '+found);
}
function test()
{

    const tests = [
    {pattern:'(a|b)*', nfa:NFA(NFA([NFA('a'), NFA('b')], '|'), '*')},
    {pattern:'^(a|b)*$', nfa:NFA([NFA('', '^'), NFA(NFA([NFA('a'), NFA('b')], '|'), '*'), NFA('', '$')], ',')},
    {pattern:'aa(b+)', nfa:NFA([NFA('aa'), NFA(NFA('b'), '+')], ',')},
    {pattern:'(a+)(b+)', nfa:NFA([NFA(NFA('a'), '+'), NFA(NFA('b'), '+')], ',')},
    {pattern:'bababa', nfa:NFA('bababa', {errors:1})},
    {pattern:'bababa', nfa:NFA('bababa', {errors:1,transpositions:true})},
    {pattern:'^(a+)(b+)$', nfa:NFA(NFA([NFA('', '^'), NFA(NFA('a'), '+'), NFA(NFA('b'), '+'), NFA('', '$')], ','), {total_errors:2})}
    ];

    let test = tests[6];
    test_case(test.nfa, test.pattern, "aaabbbbb"); // 0 errors
    test_case(test.nfa, test.pattern, "ababbbbb"); // 1 errors
    test_case(test.nfa, test.pattern, "abababbb"); // 2 errors
    test_case(test.nfa, test.pattern, "abababab"); // 3 errors
    test_case(NFA('ab'), 'ab', "ab");
    test = tests[2];
    test_case(test.nfa, test.pattern, "ababbbbaab");
    test_case(test.nfa, test.pattern, "aaababbbba");
    test = tests[3];
    test_case(test.nfa, test.pattern, "aaaaaaabbbbb");
    test = tests[4];
    test_case(test.nfa, test.pattern, "baabba");
    test = tests[5];
    test_case(test.nfa, test.pattern, "baabba");
    test = tests[4];
    test_case(test.nfa, test.pattern, "baaabbaaba");
    test = tests[5];
    test_case(test.nfa, test.pattern, "baaabbaaba");
    test = tests[4];
    test_case(test.nfa, test.pattern, "babab");
    test = tests[5];
    test_case(test.nfa, test.pattern, "babab");
    test = tests[4];
    test_case(test.nfa, test.pattern, "ababa");
    test = tests[5];
    test_case(test.nfa, test.pattern, "ababa");
    for (let i=0; i<10; ++i)
    {
        let string = create_string(['a', 'b'], 10);

        echo();
        tests.forEach(test => test_case(test.nfa, test.pattern, string));
    }
}

test();
