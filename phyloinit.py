from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="fishtank2")
db = server.new_database("fish", description="Just for testing")
server.commit()
