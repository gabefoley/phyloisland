import servers
import models
import BLAST
import utilities
import os
import time
from flask import flash
import subprocess
import Bio
from Bio import SeqIO

def getFeatureLocation(ids, reference, recordName, recordLocation):

    print (reference)

    dbpath = "files/temp_blastfiles.fasta"
    querypath = "files/psi.xml"
    reference_path = "files/blast_reference"

    SeqIO.write(reference, reference_path, "fasta")

    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
    for record in query.all():
        seq_record = servers.bio_db.lookup(primary_id=record.name)
        BLAST.makeBlastDB(seq_record, dbpath)

        while not os.path.exists(dbpath):
            time.sleep(1)

        if os.path.isfile(dbpath):
            BLAST.tBlastN(dbpath, reference_path, 0.00005)

            while not os.path.exists(querypath):
                time.sleep(1)

            if os.path.isfile(querypath):
                # with open(querypath, "w") as handle:
                #     handle.replace("b'", "")
                #     handle.replace("'","")


                blast_info = BLAST.getBlastInfo(open(querypath))

                if ("sequence" in blast_info and "location" in blast_info):

                    # Update the record with the new location
                    setattr(record, recordName, blast_info["sequence"].replace("-", ""))
                    setattr(record, recordLocation, blast_info["location"])

                    # Commit the changed record
                    servers.db.session.add(record)
                    servers.db.session.commit()



                else:
                    flash("Couldn't find an %s region in %s" % (recordName, record.name))
            # Remove created files
            utilities.removeFile(dbpath, querypath)

        else:
            raise ValueError("%s isn't a file!" % "files/temp_blastfiles.fasta")


def get_feature_location_with_profile(ids, reference, recordName, recordLocation):
    """
    Annotate a genome sequence with a feature location based on a profile
    :param ids: Genome sequences to annotate
    :param reference: Profile to annotate based on
    :param recordName: Which feature field to update
    :param recordLocation: Which feature location field to update
    :return:
    """

    dbpath = "files/temp_blastfiles.fasta"
    cleaned_path = "files/cleanedblast.fasta"

    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
    for record in query.all():
        seq_record = servers.bio_db.lookup(primary_id=record.name)
        # print ('here is seq record')
        # print (seq_record)
        # print(help(seq_record.seq.translate()))
        # print (sew)
        # print(seq_record.seq)
        # record = servers.db.lookup(accession=record.name)

        # BLAST.makeBlastDB(seq_record, dbpath)

        # pops = seq_record.seq.replace("b'", "").replace("'", "")
        nuc_seq = Bio.Seq.Seq(str(seq_record.seq).replace("b'", "").replace("'", ""))

        with open(cleaned_path, 'w') as handle:
            handle.write(">" + seq_record.name + " " + seq_record.description + "\n" + str(nuc_seq.translate(stop_symbol="M")))

            # handle.write(">" + seq_record.name + " " + seq_record.description + "\n" + str(seq_record.seq).replace("b'", "").replace("'", "").translate())
            # Bio.SeqIO.write(seq_record, handle, 'fasta')

        #Temporary measure to remove byte characters from BlastDB


        # while not os.path.exists(dbpath):
        #     time.sleep(1)
        #
        # if os.path.isfile(dbpath):
        #     textData = None
        #     with open(dbpath, "r") as handle:
        #         textData = handle.read()
        #     textData.replace("b'", "")
        #     textData.replace("'", "")
        #
        #     with open(dbpath, 'w') as handle:
        #         handle.write(cleaned_path)
        #
        while not os.path.exists(cleaned_path):
            time.sleep(1)

        if os.path.isfile(cleaned_path):
            print (reference)
            print (cleaned_path)
            result = subprocess.call(["hmmsearch", '--domtblout', 'files/output.txt', reference, cleaned_path], stdout=subprocess.PIPE)
            print ('yepr')
            print (result)

        # utilities.removeFile(cleaned_path)


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
    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
    for record in query.all():
        setattr(record, recordName, "")
        setattr(record, recordLocation, "")
        setattr(record, "Overlap", "")

        # Update the database
        servers.db.session.add(record)
        servers.db.session.commit()
