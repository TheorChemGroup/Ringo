import re

def get_matches(regex, string):
    return re.findall(regex, string)

def py_print(string):
    print(string) # ...
