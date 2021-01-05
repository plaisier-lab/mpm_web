import mysql.connector

def dbconn():
	'''
	'user': app.config['USER'],
	'password': app.config['PASS'],
	'host': app.config['HOST'],
	'port': app.config['port'],
	'database': app.config['DB'],
	'''

	'''
	'user': "root",
	'password': "root",
	'host': "db",
	'port': "3306",
	'database': "mpm_web",
	'''

	config = {
		'user': "root",
		'password': "root",
		'host': "db",
		'port': "3306",
		'database': "mpm_web",
	}
	return mysql.connector.connect(**config)