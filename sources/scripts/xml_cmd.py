#!/usr/bin/python
from sys import argv
from xml.etree.ElementTree import parse
import os

# Test if file argv[1] exists

# Open the file argv[1] and parse it
tree = parse(argv[1])
root = tree.getroot()

cmds = [];

def baseName():
	return "./build/";
#end

def libName(name):
	if(os.name == "posix"):
		return "lib" + name + ".so";
	elif(os.name == "nt"):
		return name + ".dll";
	#endif
#end

def parseFit(xmlNode):

	cmd = "";

	# Parse the command line arguments
	for child in xmlNode:
		
		if(child.tag == "function"):
			cmd += " --func " + baseName() + libName(child.attrib["plugin"]);
			cmd += " --output " + child.attrib["file"];

		elif(child.tag == "data"):
			cmd += " --input " + child.attrib["file"];
			#cmd += " --data " + baseName() + libName(child.attrib["plugin"]);

		elif(child.tag == "fitter"):
			cmd += " --fitter " + baseName() + libName(child.attrib["plugin"]);
		#endif

	#end

	return cmd;
#end

# Command lines creation
for child in root:

	cmd = "";

	if(child.tag == "fit"):
		cmd += "./build/data2brdf";
		cmd += parseFit(child);
	#end

	os.system(cmd);


#end
