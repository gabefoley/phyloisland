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
import phyloisland
import resultread
import ToxinGraphicsMain
import glob

def getFeatureLocation(ids, reference, query_name, query_location, query_length, closest_to=-1):


    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))


    for record in query.all():
        print("Looking for an %s region in %s" % (query_name, record.name))
        dbpath = "tmp/temp_blastfiles" + str(record.name).replace(" (", "_").replace(" ", "_").replace(")", "_") + ".fasta"
        reference_path = "tmp/blast_reference"
        output_path = "tmp/" + str(record.name).replace(" (", "_").replace(" ", "_").replace(")", "_") + query_name + phyloisland.randstring(5) + "_blast_results.xml"


        # Write out the reference sequnece
        SeqIO.write(reference, reference_path, "fasta")

        seq_record = servers.bio_db.lookup(primary_id=record.name)


        # Make a BLAST database of the genome to search
        BLAST.makeBlastDB(seq_record, dbpath)

        while not os.path.exists(dbpath):
            time.sleep(1)

        if os.path.isfile(dbpath):

            # Perform the BLAST search
            BLAST.tBlastN(dbpath, reference_path, output_path, 0.00005)

            while not os.path.exists(output_path):
                time.sleep(1)

            if os.path.isfile(output_path):
                print("Results of BLAST search have been written to %s \n" % output_path)

                blast_info = BLAST.getBlastInfo(open(output_path), closest_to)


                if ("sequence" in blast_info and "location" in blast_info):

                    # Update the record with the new location
                    setattr(record, query_name, blast_info["sequence"].replace("-", ""))
                    setattr(record, query_location, blast_info["location"])
                    setattr(record, query_length, len(blast_info["sequence"].replace("-", "")))

                    # Commit the changed record
                    servers.db.session.add(record)
                    servers.db.session.commit()



                else:
                    flash("Couldn't find an %s region in %s" % (query_name, record.name))

            # Remove created files
            utilities.removeFile(dbpath, dbpath + ".nhr", dbpath + ".nin", dbpath + ".nsq",  reference_path)

        else:
            raise ValueError("%s isn't a file!" % "files/temp_blastfiles.fasta")


def get_feature_location_with_profile(ids, output_dir, profile_name, recordName, recordLocation, region):
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
        id_name = record.name
        species = seq_record.annotations.get('organism')

        identifier = id_name + "_" + species
        identifier = identifier.replace(" ", "_").replace(".","_")

        # Create a path to write the translated genomic sequence to
        random_id = phyloisland.randstring(5)


        # Get the nucleotide sequence of the genome
        nuc_seq = Bio.Seq.Seq(str(seq_record.seq).replace("b'", "").replace("'", ""))
        outpath = output_dir + "/" + identifier + "/" + region +"/" + profile_name + "/"
        outpath.replace(" ", "_")
        # Check three forward reading frames
        if not os.path.exists(outpath):
            os.makedirs(outpath.replace(" ", "_"))
        
        for forward in [True, False]:
            for i in list(range(0, 3)):

                strand = "_forward_" +str(i) if forward else "_backward_" + str(i)
                sequence = nuc_seq[i:] if forward else nuc_seq.reverse_complement()[i:]

                cleaned_path = outpath + identifier + "_" + random_id + strand + "_translated_genome.fasta"
                hmmsearch_results = outpath + identifier + "_" + random_id + strand + "_hmmsearch_results.fasta"

                domScore = 100

                cleaned_path = cleaned_path.replace(" ", "_")
                hmmsearch_results = hmmsearch_results.replace(" ", "_")

                # Translate the nucleotide genome sequence to a protein sequence
                with open(cleaned_path, 'w') as handle:

                    if seq_record.name == "<unknown name>":
                        handle.write(">" + seq_record.id + " " + seq_record.annotations.get('organism') + "\n" + str(
                            sequence.translate(stop_symbol="M")))
                    else:

                        handle.write(">" + seq_record.name + " " + seq_record.description + "\n" + str(sequence.translate(stop_symbol="M")))

                print ("Writing the %s sequence with the species %s to %s" % (seq_record.id, seq_record.annotations.get('organism'), cleaned_path))

                while not os.path.exists(cleaned_path):
                    time.sleep(1)

                if os.path.isfile(cleaned_path):

                    stdoutdata = subprocess.getoutput("hmmsearch -o %s --domT %s %s %s" % (hmmsearch_results, domScore, 'tmp/' +region+"_profile.hmm", cleaned_path))

                    print (stdoutdata)

                    print ("The results from the HMM search have been written to %s \n" % hmmsearch_results)

        print("Creating a diagram of %s region hits" % (region))
        # for regions in hmmer_outputs/organism add reg to all_reg
        all_reg = []
        for infile in glob.glob(os.path.join("/" +species + '/*')):
            if '.' not in infile:
                all_reg.append(infile)
        hmmerout = []

        # add handler to HMMread for output paths 
        for reg in all_reg:
            hmmerout.append(resultread.HMMread(reg, record))
        ToxinGraphicsMain.writeHMMToImage(hmmerout, output_dir + "/" + identifier, identifier, seq_record, species)
        print("Diagram has been written to %s directory" %(output_dir))
        print("Creating a Sequence file containing all %s region hits" % (region))
	# TODO Add better invariant for Shotgun Naming
        ToxinGraphicsMain.writeHmmToSeq(hmmerout, output_dir +"/" + identifier, identifier, seq_record, species)
        print("WIP Seq may not work")


# Function to be called when wanting to generate Genome Diagram and GenBank output
def generateOutput(ids, output_dir, diagram=True, genbank=True, fasta=False, expand=False):

    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
    fasta_region_dict = {}
    for record in query.all():
        seq_record = servers.bio_db.lookup(primary_id=record.name)
        id_name = record.name
        species = seq_record.annotations.get('organism')

        identifier = id_name + "_" + species
        identifier = identifier.replace(" ", "_").replace(".","_")

        all_reg = []

        for infile in glob.glob(os.path.join(output_dir + "/" + identifier + '/*')):
            if '.' not in infile:
                all_reg.append(infile)

        outpath = os.path.join(output_dir + "/" + id_name.replace(" ", "_") + "_" + species.replace(" ", "_"))

        hmmerout = []

        for reg in all_reg:

            hmmerout.append(resultread.HMMread(reg, record, expand))

            if fasta:
                reg_name = reg.split("/")[-1]

                if reg_name not in fasta_region_dict:
                    fasta_region_dict[reg_name] = []


                # For a given region, create the sequence files
                fasta_output = utilities.createFASTAFromHMMOutput(seq_record, hmmerout, record.name, reg_name)

                # Write the region sequence files out to a genome specific FASTA file
                fasta_outpath = "%s/%s%s.fasta" % (outpath, reg_name, "_expanded" if expand else "")
                utilities.saveFASTA(fasta_output, fasta_outpath)

                # Add the region to a dictionary so we can write out a non-genome specific FASTA file
                fasta_region_dict[reg_name] = fasta_region_dict[reg_name] + fasta_output

                hmmerout = []

        if diagram:
            print ('here')
            ToxinGraphicsMain.writeHMMToImage(hmmerout, outpath, identifier, seq_record,
                                          species, expand)
            print("Diagram has been written to %s directory" % (output_dir))

        if genbank:
            # TODO Add better invariant for Shotgun Naming
            ToxinGraphicsMain.writeHmmToSeq(hmmerout, outpath, identifier, seq_record, species, expand)

    # Write out the non-genome specific FASTA file/s
    if fasta:
        for reg_name, output in fasta_region_dict.items():
            outpath = "%s/hmm_outputs/%s_%s%s.fasta" % (os.getcwd(), reg_name, "expanded_" if expand else "", str(len(ids)))
            utilities.saveFASTA(output, outpath)






def deleteFeature(**kwargs):
    """
    Delete a certain feature in database
    :param ids: List of ids to delete features from
    :param recordName: Name of feature to delete
    :param recordLocation: Name of feature location to delete
    :return:
    """
    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(kwargs["ids"]))
    for record in query.all():

        for key, value in kwargs.items():
            if key == "string":
                for item in value:
                    setattr(record, item, "")
            if key == "int":
                for item in value:
                    setattr(record, item, None)

            elif key == "bool":
                for item in value:
                    setattr(record, item, 0)

        # Update the database
        servers.db.session.add(record)
        servers.db.session.commit()

