from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="fishtank_stable")
db = server.new_database("stable", description="Just for testing")
server.commit()
