"use strict";

const Matchy = require('../../src/js/Matchy.js');
const echo = console.log;

function create_string(alphabet, n)
{
    let s = '';
    for (let i=0; i<n; ++i)
    {
        s += alphabet[Math.round(Math.random()*(alphabet.length-1))];
    }
    return s;
}
function create_pattern(string, n)
{
    const i = Math.round(Math.random()*(string.length-n));
    return string.slice(i, i+n);
}
function test_case(matchy, algorithm, pattern, string, offset)
{
    offset = offset || 0;
    const matcher = matchy[algorithm](pattern);
    const found = matcher(string, offset);
    const index = string.indexOf(pattern, offset);
    echo(algorithm+'("'+pattern+'", "'+string+'", '+offset+') = '+found+(found === index ? ' (true)' : ' (expected '+index+')'));
}
function test()
{
    const matchy = new Matchy();

    const algorithms = [
    'fsa',
    'rabinkarp',
    'knuthmorrispratt',
    'twoway',
    'boyermoore'
    ];

    for (let i=0; i<10; ++i)
    {
        let string = create_string(['a', 'b'/*, 'c', 'd'*/], 10);
        let pattern = create_pattern(string, 5);

        echo();
        algorithms.forEach(algorithm => test_case(matchy, algorithm, pattern, string));
    }
    /*
    // problematic
    echo();
    echo('problematic');
    ([
    ['twoway', "babab", "aababababb"],
    ['twoway', "babab", "baaabababb"],
    ['boyermoore', "babab", "aababababb"],
    ['boyermoore', "abaaa", "baabaaabab"]
    ]).forEach(entry => test_case(matchy, entry[0], entry[1], entry[2]));
    */
}

test();
