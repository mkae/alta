#!/usr/bin/python

"""@package Python command creator package

With this package, you can generate ALTA command lines within python.
"""

import os


class command:

	def __init__(self):
		self.name = 'data2brdf'
	#end

	def run(self):

		cmd = self.name
		for opt in self.options:
			cmd + ' ' + opt
		#endfor

		os.system(cmd)
	#end

	def set_input(self, filename):
		self.options.prepend('--input ' + filename)
	#end
	
	def set_output(self, filename):
		self.options.prepend('--output ' + filename)
	#end

	name    = ''
	options = []
#endclass
