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

def getFeatureLocation(ids, reference, record_name, record_location, record_length):

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
                    setattr(record, record_name, blast_info["sequence"].replace("-", ""))
                    setattr(record, record_location, blast_info["location"])
                    setattr(record, record_length, len(blast_info["sequence"].replace("-", "")))

                    # Commit the changed record
                    servers.db.session.add(record)
                    servers.db.session.commit()



                else:
                    flash("Couldn't find an %s region in %s" % (record_name, record.name))
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

    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
    for record in query.all():
        seq_record = servers.bio_db.lookup(primary_id=record.name)

        # Create a path to write the translated genomic sequence to
        cleaned_path = "tmp/" + seq_record.id + "_translated_genome.fasta"
        hmmsearch_results = "tmp/" + seq_record.id + "_hmmsearch_results.fasta"

        # Get the nucleotide sequence of the genome
        nuc_seq = Bio.Seq.Seq(str(seq_record.seq).replace("b'", "").replace("'", ""))

        # Translate the nucleotide genome sequence to a protein sequence
        with open(cleaned_path, 'w') as handle:
            handle.write(">" + seq_record.name + " " + seq_record.description + "\n" + str(nuc_seq.translate(stop_symbol="M")))

        print ("Writing the %s sequence with the species %s to %s \n" % (seq_record.id, seq_record.annotations.get('organism'), cleaned_path))

        while not os.path.exists(cleaned_path):
            time.sleep(1)

        if os.path.isfile(cleaned_path):

            stdoutdata = subprocess.getoutput("hmmsearch -o %s %s %s" % (hmmsearch_results, reference, cleaned_path))
            print ("The results from the HMM search have been written to %s" % hmmsearch_results)
            # result = subprocess.call(["hmmsearch", 'files/output.txt', reference, cleaned_path], stdout=subprocess.PIPE)
            # for x in result:
            #     print (x)

        # utilities.removeFile(reference, cleaned_path)


def checkFeatureAlreadyAnnotated():
    return

def annotateNewFeature():
    return

def addToDatabase():
    return




def deleteFeature(ids, record_name, record_location, record_length):
    """
    Delete a certain feature in database
    :param ids: List of ids to delete features from
    :param recordName: Name of feature to delete
    :param recordLocation: Name of feature location to delete
    :return:
    """
    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
    for record in query.all():
        setattr(record, record_name, "")
        setattr(record, record_location, "")
        setattr(record, record_length, None)
        setattr(record, "overlap", "")
        setattr(record, "distance", "")


        # Update the database
        servers.db.session.add(record)
        servers.db.session.commit()
