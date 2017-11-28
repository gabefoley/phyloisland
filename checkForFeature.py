import servers
import models
import BLAST
import utilities
import os
import time
from flask import flash

def getFeatureLocation(ids, reference, recordName, recordLocation):

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
                blast_info = BLAST.getBlastInfo(open(querypath))

                if ("sequence" in blast_info and "location" in blast_info):

                    # Update the record with the new location
                    setattr(record, recordName, blast_info["sequence"])
                    setattr(record, recordLocation, blast_info["location"])

                    # Commit the changed record
                    servers.db.session.add(record)
                    servers.db.session.commit()


                    # Remove created files
                    utilities.removeFile(dbpath, querypath)
                else:
                    flash("Couldn't find an %s region in %s" % (recordName, record.name))

        else:
            raise ValueError("%s isn't a file!" % "files/temp_blastfiles.fasta")


def checkFeatureAlreadyAnnotated():
    return

def annotateNewFeature():
    return

def addToDatabase():
    return




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




