#!/usr/local/bin/python
__author__ = 'miceli'

file = open("errors.txt", 'r', -1)

#78
i = 0
for line in file.readlines():
    print line

