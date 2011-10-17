#py grep command
#sample command: grep("^x",dir())
#syntax: grep(regexp_string,list_of_strings_to_search)

import re

def grep(string,list):
    expr = re.compile(string)
    for text in list:
        match = expr.search(text)
        if match != None:
            print match.string

# Two alternate solutions, courtesy of Gene H
# These are more elegant, and more detailed descriptions can be found here:
# http://www.faqs.org/docs/diveintopython/apihelper_filter.html
# http://www.faqs.org/docs/diveintopython/regression_filter.html

def grep(string,list):
    expr = re.compile(string)
    return [elem for elem in list if expr.match(elem)]

def grep(string,list):
    expr = re.compile(string)
    return filter(expr.search,list)

# Note a subtle difference between the above two:
# because the first uses expr.match(), it will only return exact matches
# while the latter will return a list containing strings that contain the
# search string in any position
# e.g.:
# list = ['normalize','size','nonzero','zenith']
# grep('ze',list)
# First one returns: ['zenith']
# Second one returns: ['normalize','size','nonzero','zenith']

def grepv(string,list):
    """ grep -v - return elements that do NOT contain the string """
    expr = re.compile(string)
    return [elem for elem in list if not expr.search(elem)]
