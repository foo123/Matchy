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

def test_case(nfa, pattern, string, match, errors):
    found = nfa.match(string, 0, False, True)
    print('fuzzynfa("'+pattern+'", "'+string+'") = '+str(found['match'])+', errors '+str(found['errors'])+' (expected '+str(match)+', errors '+str(errors)+')')

def test():
    test = {'pattern':'^(a+)(b+)$', 'nfa':NFA(NFA([NFA('', '^'), NFA(NFA('a'), '+'), NFA(NFA('b'), '+'), NFA('', '$')], ','), {'total_errors':2,'word_level':False})}
    test_case(test['nfa'], test['pattern'], "aaabbbbb", 0, 0) # 0 errors
    test_case(test['nfa'], test['pattern'], "ababbbbb", 0, 1) # 1 errors
    test_case(test['nfa'], test['pattern'], "abababbb", 0, 2) # 2 errors
    test_case(test['nfa'], test['pattern'], "abababab", -1, 3) # 3 errors
    test_case(test['nfa'], test['pattern'], "aabababbbb", 0, 2) # 2 errors
    test_case(test['nfa'], test['pattern'], "baabaaabbb", 0, 2) # 2 errors
    print()

    test = {'pattern':'(aaa)(bbb)', 'nfa':NFA(NFA([NFA('aaa'), NFA('bbb')], ','), {'total_errors':1,'word_level':True})}
    test_case(test['nfa'], test['pattern'], "aaabbb", 0, 0) # 0 errors
    test_case(test['nfa'], test['pattern'], "bbbaaa", 0, 1) # 1 errors, deletion
    test_case(test['nfa'], test['pattern'], "bbb", 0, 1) # 1 errors, deletion
    test_case(test['nfa'], test['pattern'], "cbbb", 0, 1) # 1 errors, deletion or substitution
    test_case(test['nfa'], test['pattern'], "aacbbb", 0, 1) # 1 errors, deletion or substitution
    test_case(test['nfa'], test['pattern'], "aaacbbb", 0, 1) # 1 errors, insertion
    print()

    test = {'pattern':'^(aaa)(bbb)$', 'nfa':NFA(NFA([NFA('', '^'), NFA('aaa'), NFA('bbb'), NFA('', '$')], ','), {'total_errors':1,'word_level':True,'transpositions':True})}
    test_case(test['nfa'], test['pattern'], "aaabbb", 0, 0) # 0 errors
    test_case(test['nfa'], test['pattern'], "bbbaaa", 0, 1) # 1 errors, transposition
    test_case(test['nfa'], test['pattern'], "bbb", 0, 1) # 1 errors, deletion
    test_case(test['nfa'], test['pattern'], "cbbb", 0, 1) # 1 errors, substitution
    test_case(test['nfa'], test['pattern'], "aacbbb", 0, 1) # 1 errors, substitution
    test_case(test['nfa'], test['pattern'], "aaacbbb", 0, 1) # 1 errors, insertion
    test_case(test['nfa'], test['pattern'], "cccbbb", 0, 1) # 1 errors, substitution
    print()

test()