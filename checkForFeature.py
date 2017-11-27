
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO, pairwise2
import servers
import models
from flask import flash
import re
import BLAST
import utilities
import os
import time

def getFeatureLocation(ids, reference):
    query = models.SequenceRecords.query.filter(models.SequenceRecords.uid.in_(ids))
    for record in query.all():
        seq_record = servers.bio_db.lookup(primary_id=record.name)
        BLAST.makeBlastDB(seq_record)

        while not os.path.exists("files/temp_blastfiles.fasta"):
            time.sleep(1)

        if os.path.isfile("files/temp_blastfiles.fasta"):
            BLAST.tBlastN("files/temp_blastfiles.fasta", reference, 0.00005)

            while not os.path.exists("files/psi.xml"):
                time.sleep(1)

            if os.path.isfile("files/psi.xml"):
                print ('made')

                location = BLAST.getBlastLocation(open("files/psi.xml"))
                return location

                utilities.removeFile("files/temp_blastfiles.fasta", "files/psi.xml")
        else:
            raise ValueError("%s isn't a file!" % "files/temp_blastfiles.fasta")


def checkFeatureAlreadyAnnotated():
    return

def annotateNewFeature():
    return

def addToDatabase():
    return




# def checkFeature(ids, reference, recordName, recordLocation):
#
#     query = models.SequenceRecords.query.filter(models.SequenceRecords.uid.in_(ids))
#     for record in query.all():
#         best_align = 0
#         best_location = None
#         best_seq = None
#         print(record.name)
#         seq_record = servers.bio_db.lookup(primary_id=record.name)
#
#         # print(seq_record.description)
#         # print(seq_record)
#
#         for feature in seq_record.features:
#             # print(feature)
#             if 'translation' in feature.qualifiers:
#
#                 alignment = pairwise2.align.globalms(reference, feature.qualifiers['translation'][0],
#                                                      2, -1, -1, -.5, score_only=True)
#                 if alignment > best_align:
#                     best_align = alignment
#                     best_seq = feature.qualifiers['translation'][0]
#                     best_location = feature.location
#
#         if (best_seq):
#             print("best align is ", best_align)
#             flash("Found an A1 region in %s" % record.name)
#             location = re.search(r"\d*:\d*", str(best_location))
#
#             setattr(record, recordName, best_seq)
#             print('location is')
#             print(location)
#             if location:
#                 print(location.group(0))
#                 setattr(record, recordLocation, location.group(0))
#                 servers.db.session.add(record)
#                 servers.db.session.commit()
#
#         # alignment = pairwise2.align.globalms(yenA1, best_seq,
#         #                                      2, -1, -1, -.5)


def deleteFeature(ids, recordName, recordLocation):
    query = models.SequenceRecords.query.filter(models.SequenceRecords.uid.in_(ids))
    for record in query.all():
        setattr(record, recordName, "")
        setattr(record, recordLocation, "")
        setattr(record, "Overlap", "")
        servers.db.session.add(record)
        servers.db.session.commit()




