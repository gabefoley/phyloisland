import servers
import models
import BLAST
import utilities
import os
import time

def getFeatureLocation(ids, reference):

    dbpath = "files/temp_blastfiles.fasta"
    querypath = "files/psi.xml"

    query = models.SequenceRecords.query.filter(models.SequenceRecords.uid.in_(ids))
    for record in query.all():
        print (record)
        seq_record = servers.bio_db.lookup(primary_id=record.name)
        BLAST.makeBlastDB(seq_record, dbpath)

        while not os.path.exists(dbpath):
            time.sleep(1)

        if os.path.isfile(dbpath):
            BLAST.tBlastN(dbpath, reference, 0.00005)

            while not os.path.exists(querypath):
                time.sleep(1)

            if os.path.isfile(querypath):
                location = BLAST.getBlastLocation(open(querypath))
                print (location )
                # return location

                utilities.removeFile(dbpath, querypath)
                return location
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
    """
    Delete a certain feature in database
    :param ids: List of ids to delete features from
    :param recordName: Name of feature to delete
    :param recordLocation: Name of feature location to delete
    :return:
    """
    query = models.SequenceRecords.query.filter(models.SequenceRecords.uid.in_(ids))
    for record in query.all():
        setattr(record, recordName, "")
        setattr(record, recordLocation, "")
        setattr(record, "Overlap", "")

        # Update the database
        servers.db.session.add(record)
        servers.db.session.commit()




