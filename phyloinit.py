from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="bioseqdb")
db = server.new_database("grass", description="Just for testing")
server.commit()
