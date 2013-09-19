#!/usr/bin/python
from sys import argv
from xml.etree.ElementTree import parse

# Test if file argv[1] exists

# Open the file argv[1] and parse it
tree = parse(argv[1])
root = tree.getroot()

# Command lines creation
for child in root:
	print child.tag, child.attrib
