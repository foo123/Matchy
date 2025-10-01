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

def create_string(alphabet, n):
    s = '';
    for i in range(n):
        s += alphabet[random.randint(0, len(alphabet)-1)];
    return s

def create_pattern(string, n):
    i = random.randint(0, len(string)-n)
    return string[i:i+n]

def test_case(matchy, algorithm, pattern, string, offset = 0):
    matcher = getattr(matchy, algorithm)(pattern)
    found = matcher(string, offset)
    index = string.find(pattern, offset)
    print(algorithm+'("'+pattern+'", "'+string+'", '+str(offset)+') = '+str(found)+(' (true)' if found == index else ' (expected '+str(index)+')'))

def test():
    matchy = Matchy()

    algorithms = [
    'fsa',
    'rabinkarp',
    'knuthmorrispratt',
    'twoway',
    'boyermoore'
    ]

    for i in range(10):
        string = create_string(['a', 'b'], 10) #['a', 'b', 'c', 'd']
        pattern = create_pattern(string, 5)

        print()
        for algorithm in algorithms:
            test_case(matchy, algorithm, pattern, string)

    ## problematic
    #print()
    #print('problematic')
    #problematic = [
    #['twoway', "babab", "aababababb"],
    #['twoway', "babab", "baaabababb"],
    #['boyermoore', "babab", "aababababb"],
    #['boyermoore', "abaaa", "baabaaabab"]
    #]
    #for entry in problematic: test_case(matchy, entry[0], entry[1], entry[2])

test()