# -*- coding: utf-8 -*-
import os, sys

DIR = os.path.dirname(os.path.abspath(__file__))

def import_module(name, path):
    #import imp
    #try:
    #    mod_fp, mod_path, mod_desc  = imp.find_module(name, [path])
    #    mod = getattr( imp.load_module(name, mod_fp, mod_path, mod_desc), name )
    #except ImportError as exc:
    #    mod = None
    #    sys.stderr.write("Error: failed to import module ({})".format(exc))
    #finally:
    #    if mod_fp: mod_fp.close()
    #return mod
    import importlib.util, sys
    spec = importlib.util.spec_from_file_location(name, path+name+'.py')
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return getattr(mod, name)

# import the Matchy.py (as a) module, probably you will want to place this in another dir/package
Matchy = import_module('Matchy', os.path.join(DIR, '../../src/py/'))
if not Matchy:
    print ('Could not load the Matchy Module')
    sys.exit(1)

import random

def NFA(input, type = 'l'):
    return Matchy.NFA(input, type)

def create_string(alphabet, n):
    s = '';
    for i in range(n):
        s += alphabet[random.randint(0, len(alphabet)-1)];
    return s

def test_case(nfa, pattern, string, offset = 0):
    found = nfa.match(string, offset)
    print('nfa("'+pattern+'", "'+string+'", '+str(offset)+') = '+str(found[0])+', errors '+str(found[1]))

def test():

    tests = [
    {'pattern':'(a|b)*', 'nfa':NFA(NFA([NFA('a'), NFA('b')], '|'), '*')},
    {'pattern':'^(a|b)*$', 'nfa':NFA([NFA('', '^'), NFA(NFA([NFA('a'), NFA('b')], '|'), '*'), NFA('', '$')], ',')},
    {'pattern':'aa(b+)', 'nfa':NFA([NFA('aa'), NFA(NFA('b'), '+')], ',')},
    {'pattern':'(a+)(b+)', 'nfa':NFA([NFA(NFA('a'), '+'), NFA(NFA('b'), '+')], ',')},
    {'pattern':'bababa', 'nfa':NFA('bababa', {'errors':1})},
    {'pattern':'bababa', 'nfa':NFA('bababa', {'errors':1,'transpositions':True})}
    ]

    for i in range(10):
        string = create_string(['a', 'b'], 10)

        print()
        for test in tests:
            test_case(test['nfa'], test['pattern'], string)

    print()
    # some special test cases
    test_case(NFA('ab'), 'ab', "ab")
    test = tests[2]
    test_case(test['nfa'], test['pattern'], "ababbbbaab")
    test_case(test['nfa'], test['pattern'], "aaababbbba")
    test = tests[3]
    test_case(test['nfa'], test['pattern'], "aaaaaaabbbbb")
    test = tests[4]
    test_case(test['nfa'], test['pattern'], "baabba")
    test = tests[5]
    test_case(test['nfa'], test['pattern'], "baabba")
    test = tests[4]
    test_case(test['nfa'], test['pattern'], "baaabbaaba")
    test = tests[5]
    test_case(test['nfa'], test['pattern'], "baaabbaaba")
    test = tests[4]
    test_case(test['nfa'], test['pattern'], "babab")
    test = tests[5]
    test_case(test['nfa'], test['pattern'], "babab")
    test = tests[4]
    test_case(test['nfa'], test['pattern'], "ababa")
    test = tests[5]
    test_case(test['nfa'], test['pattern'], "ababa")

test()