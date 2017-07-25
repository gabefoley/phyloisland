from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="hedgehogdb")
db = server.new_database("hedgehog", description="Just for testing")
server.commit()
