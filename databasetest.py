from BioSQL import BioSeqDatabase
from Bio import Entrez, SeqIO
server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="bioseqdb")


# db = server["phylotest"]
#
# for key in db:
#     record = db[key]
#     del db[key]
# server.commit()

# Add new sub-database
db = server.new_database("phylotest", description="Test environment for Phylo Island")


#
# for key, record in db.items():
#     print (key, record)


# for identifier in ['AF189967.1'] :
#     seq_record = db.lookup(gi=identifier)
#     print (seq_record.id)
#     print (seq_record.annotations["GABE"])

# print ("This database contains %i records" % len(db))
# for key, record in db.items():
#     for annot in record.annotations:
#         print (annot)
#     print ("Key %r maps to a sequence record with id %s and annotations %s" % (key, record.id, record.annotations))




# handle = Entrez.efetch(db="nuccore", id="6273244", rettype="gb", email="gabriel.foley@uqconnect.edu.au", retmode="text")
#
# record_dict = list(SeqIO.parse(handle, "gb"))
#
# for record in record_dict:
#     print (record.id)
#     if record.id == "AF189967.1":
#         record.annotations["GABE"] = "Added this"
#     for annot in record.annotations:
#         print (annot)
#
# count = db.load(record_dict)
# server.commit()

#
# count = db.load(records)
# print ("Loaded %i records" % count)
server.commit()