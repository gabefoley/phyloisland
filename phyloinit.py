from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="walrus")
db = server.new_database("water", description="Just for testing")
server.commit()
