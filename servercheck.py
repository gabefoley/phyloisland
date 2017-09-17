from servers import *
name = "name"
species = "species"
strain = "strain"
description = "description"
description = "desc"
a1 = "Not tested"
a1_loc =  "Not tested"
a2 =  "Not tested"
a2_loc =  "Not tested"
overlap =  "Not tested"
distance = "Not tested"
sequence = "AAPPG"



class SequenceRecords(db.Model):
    __tablename__ = 'seqrecord'
    uid = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text)
    species = db.Column(db.String(255))
    strain = db.Column(db.String(255))
    description = db.Column(db.Text)
    a1 = db.Column(db.Text)
    a1_loc = db.Column(db.Text)
    a2 = db.Column(db.Text)
    a2_loc = db.Column(db.Text)

    overlap = db.Column(db.String(255))
    distance = db.Column(db.VARCHAR(255))
    sequence = db.Column(db.Text)

    def __init__(self, name="", species="", strain="", description="", a1="", a1_loc="", a2="", a2_loc="", overlap="",
                 distance="", sequence=""):
        self.name = name
        self.species = species
        self.strain = strain
        self.description = description
        self.a1 = a1
        self.a1_loc = a1_loc
        self.a2 = a2
        self.a2_loc = a2_loc
        self.overlap = overlap
        self.distance = distance
        self.sequence = sequence

entry = SequenceRecords(name, species, strain, description, a1, a1_loc, a2, a2_loc, overlap, distance, sequence)

db.session.add(entry)
db.session.commit()
# bio_server.commit()