from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="bioseqdb2")
db = server.new_database("phylomain", description="Just for testing")
server.commit()
