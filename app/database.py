import mysql.connector
import os

def dbconn():
	host = "db"
	if os.environ.get("MPM_WEB_SQL_HOST") != None:
		host = os.environ.get("MPM_WEB_SQL_HOST")
	
	port = "3306"
	if os.environ.get("MPM_WEB_SQL_PORT") != None:
		port = os.environ.get("MPM_WEB_SQL_PORT")
	
	config = {
		'user': "root",
		'password': "root",
		'host': host,
		'port': port,
		'database': "mpm_web",
	}
	return mysql.connector.connect(**config)