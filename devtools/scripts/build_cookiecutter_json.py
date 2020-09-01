import sys
import json


data = {'name': 'openforcefield',
 'channel': 'omnia',
 'python': ['3.6', '3.7'],
 'platform': ['linux-64', 'osx-64'],
 'release': sys.argv[1]}

with open('new_cookiecutter.json', 'w') as fp:
    json.dump(data, fp)
