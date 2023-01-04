import json


def read_json(filename):
    with open(filename, 'r') as file:
        dictionary = json.loads(file.read())

    return dictionary


def write_json(filename, dictionary):
    with open(filename, 'w') as file:
        file.write(json.dumps(dictionary))
