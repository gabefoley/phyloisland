from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="primates")
db = server.new_database("newsmall", description="Just for testing")
server.commit()
