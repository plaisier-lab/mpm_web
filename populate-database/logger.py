# https://stackoverflow.com/a/616672 for debugging
import sys
class Logger(object):
	def __init__(self):
		self.terminal = sys.stdout
		open("console.log", "w").close() # clear the file
		self.log = open("console.log", "a")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)

sys.stdout = Logger()