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
	return './build/';
#end

def libName(name):
	if(os.name == 'posix'):
		return 'lib' + name + '.so';
	elif(os.name == 'nt'):
		return name + ".dll";
	#endif
#end

## Relative directories
lib_dir = '';
dat_dir = '';
out_dir = '';


## Parse the configuration part of the XML file, this will set the global
## parameters such as the relative directories.
##
def parseConfiguration(xmlNode):
	global lib_dir, dat_dir, out_dir;

	for param in xmlNode.findall('parameter'):

		if(param.attrib['name'] == 'lib-dir'):
			lib_dir = param.attrib['value'];

		elif(param.attrib['name'] == 'dat-dir'):
			dat_dir = param.attrib["value"];

		elif(param.attrib['name'] == 'out-dir'):
			dat_dir = param.attrib["value"];

		#end
	#end
#end


def parseAction(xmlNode):
	global lib_dir, dat_dir, out_dir;
	cmd = '';

	inputNode = xmlNode.find('input');
	if not (inputNode is None):
		cmd += ' --input ' + inputNode.attrib['name'];
	#end
	
	outputNode = xmlNode.find('output');
	if not(outputNode is None):
		cmd += ' --output ' + outputNode.attrib['name'];
	#end

	for plugin in xmlNode.findall('plugin'):
		cmd += ' --' + plugin.attrib['type'];
		cmd += ' ' + lib_dir + '/' + libName(plugin.attrib['name']);
	#end
	
	for param in xmlNode.findall('parameter'):
		cmd += ' --' + param.attrib['name'];
		cmd += ' ' + param.attrib['value'];
	#end

	return cmd;
#end


## Parse the configuration part of the file
##
conf = root.find("configuration");
if conf is None:
	print '<<PYTHON>> no configuration specified in the XML file';
else:
	parseConfiguration(conf);
#end


## Command lines creation
##
for child in root.findall('action'):

	# Create the cmd string with the command name
	cmd = lib_dir + '/' + child.attrib['name'];

	# Parse the action
	cmd += parseAction(child);

	#print cmd;
	ret = os.system(cmd);
	if(ret != 0):
		print '<<PYTHON>> the action was not performed';
	#end


#end
