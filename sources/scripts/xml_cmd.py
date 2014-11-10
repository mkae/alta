#!/usr/bin/python


"""@package Python command creator

With this package, you can generate ALTA command lines with XML formated
files. The documentation of the XML format can be found in the Format
page.
"""

import sys
from sys import argv
from xml.etree.ElementTree import parse
import subprocess
import os
from os import listdir
from os.path import isfile, join

""" Color capable outputing
"""
BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

#following from Python cookbook, #475186
def has_colours(stream):
	if not hasattr(stream, "isatty"):
		return False
	if not stream.isatty():
		return False # auto color only on TTYs
	try:
		import curses

		curses.setupterm()
		return curses.tigetnum("colors") > 2
	except:
		# guess false in case of error
		return False


has_colours = has_colours(sys.stdout)


def printout(text, colour=WHITE):
	if has_colours:
		seq = "\x1b[1;%dm" % (30 + colour) + text + "\x1b[0m"
		sys.stdout.write(seq)
	else:
		sys.stdout.write(text)
	sys.stdout.flush()





# Test if file argv[1] exists
if len(argv) != 2:
	print '<<PYTHON>> you did not specify the input argument correctly:'
	print '           xml_cmd file.xml'
	exit()
#end


# Open the file argv[1] and parse it
tree = parse(argv[1])
root = tree.getroot()


def baseName():
	return './build/'

#end

def libName(name):
	if sys.platform == 'darwin':
		return 'lib' + name + '.dylib'
	elif os.name == 'nt':
		return name + ".dll" 
	else:
		return 'lib' + name + '.so'
	#endif

#end

##Directory delimiter
dd = '/'
if os.name == 'nt':
	dd = '\\'

## Relative directories
lib_dir = ''
dat_dir = ''
out_dir = ''


## Parse the configuration part of the XML file, this will set the global
## parameters such as the relative directories.
##
def parseConfiguration(xmlNode):
	global lib_dir, dat_dir, out_dir

	for param in xmlNode.findall('parameter'):

		if (param.attrib['name'] == 'lib-dir'):
			lib_dir = param.attrib['value']

		elif (param.attrib['name'] == 'dat-dir'):
			dat_dir = param.attrib["value"]

		elif (param.attrib['name'] == 'out-dir'):
			dat_dir = param.attrib["value"]

		#end
	#end
#end

def parseFunction(xmlNode):
	global lib_dir

	cmd = ''
	#cmd = lib_dir + dd + libName(xmlNode.attrib['name'])
	cmd = xmlNode.attrib['name']

	# Parse the parameters of the function
	for param in xmlNode.findall('parameter'):
		cmd += ' --' + param.attrib['name'] + ' ' + param.attrib['value']
	#end

	fresnel = xmlNode.find('fresnel')
	if not (fresnel is None):
		cmd += ' --fresnel ' + lib_dir + dd + libName(fresnel.attrib['name'])
	#end

	return cmd

#end

def parseFunctions(xmlNodes):
	list_len = len(xmlNodes)
	if list_len == 0:
		return ''
	elif list_len == 1:
		return ' --func ' + parseFunction(xmlNodes[0])
	else:
		cmd = ' --func ['
		for index in range(0, list_len):
			cmd += parseFunction(xmlNodes[index])
			if index != list_len - 1:
				cmd += ', '
			#end
		#end
		cmd += ']'
		return cmd
	#end

#end

def parseAction(xmlNode):
	global lib_dir, dat_dir, out_dir
	cmd = ''

	cmd += parseFunctions(xmlNode.findall('function'))

	for plugin in xmlNode.findall('plugin'):
		cmd += ' --' + plugin.attrib['type']
		#cmd += ' ' + lib_dir + dd + libName(plugin.attrib['name'])
                cmd += ' ' + plugin.attrib['name']
	#end

	for param in xmlNode.findall('parameter'):
		cmd += ' --' + param.attrib['name']
		cmd += ' ' + param.attrib['value']
	#end

        # Search the input file to the command. It is possible for the xml
        # command to loop over a list of file present in a directory using
        # the directory attrib. In that case, cmd is duplicated into a cmd
        # list.
	inputNode = xmlNode.find('input')
	outputNode = xmlNode.find('output')
	if not (inputNode is None) and not (outputNode is None):
            if 'directory' in inputNode.attrib:
                dirname = inputNode.attrib['directory']
                cmds = [cmd + ' --input ' + join(dirname, f) + ' --output ' + f + '.out' for f in listdir(dirname) if isfile(join(dirname, f))]
                cmd = cmds;
            else:
	        cmd += ' --input '  + inputNode.attrib['name']
		cmd += ' --output ' + outputNode.attrib['name']
                cmd = [cmd]
            #end
	#end

        return cmd

#end


## Parse the configuration part of the file
##
conf = root.find("configuration")
if not (conf is None):
    parseConfiguration(conf)


## Command lines creation
##
for child in root.findall('action'):

	# Create the cmd string with the command name
        cmd = ''
        if lib_dir != '':
	    cmd = lib_dir + dd + child.attrib['name']
        else:
	    cmd = child.attrib['name']


	# Parse the action
        try:
            for act in parseAction(child):

                printout(cmd + act, GREEN)
                print

	        ret = subprocess.check_call(cmd + act, shell=True)
                if ret != 0:
	            printout('<<PYTHON>> the action was not performed', RED)
	        #end
	        print '\n'
            #end
        except:
            printout('<<PYTHON>> exception caught, stopping', RED)
            print '\n'
            exit(0)
#end
