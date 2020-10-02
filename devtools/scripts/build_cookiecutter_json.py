import sys
import json


release_tag = sys.argv[1]
python_version = sys.argv[2]
ci_os = sys.argv[3]

platform_mapping = {
    "ubuntu-latest": "linux-64",
    "macOS-latest": "osx-64",
}

data = {
    "name": "openforcefield",
    "channel": "omnia",
    "python": [python_version],
    "platform": [platform_mapping[ci_os]],
    "release": release_tag,
}

with open("new_cookiecutter.json", "w") as fp:
    json.dump(data, fp)
