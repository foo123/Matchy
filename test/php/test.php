<?php
include(dirname(__FILE__) . '/../../src/php/Matchy.php');

function create_string($alphabet, $n)
{
    $s = '';
    for ($i=0; $i<$n; ++$i)
    {
        $s .= $alphabet[rand(0, count($alphabet)-1)];
    }
    return $s;
}
function create_pattern($string, $n)
{
    $i = rand(0, strlen($string)-$n);
    return substr($string, $i, $n);
}
function test_case($matchy, $algorithm, $pattern, $string, $offset = 0)
{
    $matcher = $matchy->{$algorithm}($pattern);
    $found = $matcher($string, $offset);
    $index = strpos($string, $pattern, $offset);
    if (false === $index) $index = -1;
    echo($algorithm.'("'.$pattern.'", "'.$string.'", '.$offset.') = '.$found.($found === $index ? ' (true)' : ' (expected '.$index.')')."\n");
}
function test()
{
    $matchy = new Matchy();

    $algorithms = [
    'fsa',
    'rabinkarp',
    'knuthmorrispratt',
    'twoway',
    'boyermoore'
    ];

    for ($i=0; $i<10; ++$i)
    {
        $string = create_string(['a', 'b'/*, 'c', 'd'*/], 10);
        $pattern = create_pattern($string, 5);

        echo("\n");
        foreach ($algorithms as $algorithm)
        {
            test_case($matchy, $algorithm, $pattern, $string);
        }
    }

    /*
    // problematic
    echo("\n");
    echo('problematic'."\n");
    $problematic = [
    ['boyermoore', "bbbbb", "aabaabbbbb"],
    ['twoway', "babab", "aababababb"],
    ['twoway', "babab", "baaabababb"],
    ['boyermoore', "babab", "aababababb"],
    ['boyermoore', "abaaa", "baabaaabab"]
    ];
    foreach ($problematic as $entry) test_case($matchy, $entry[0], $entry[1], $entry[2]);
    */
}

test();