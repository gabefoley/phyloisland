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
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
import phyloisland

def writeImageToFile(ids):

    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
    # Initiation of GenomeDiagram assets"
    name = "GenomeDiagram_"+phyloisland.randstring(5)
    gd_diagram = GenomeDiagram.Diagram(name)
    max_len = 0
    output_path = "tmp/"+name+".png"
    for record in query.all():
        
        """ Create a dictionary for key = feature type -> value = location """
        locs = {}
        # Prepare for literally the worst code in existence
        """ We have to pull the features from location data in database
        instead of directly from database because of current limitations """
        if getattr(record, "a1_loc") != "":
            locs["a1"] = getattr(record, "a1_loc").split(":")
        if getattr(record, "a2_loc") != "":
            locs["a2"] = getattr(record, "a2_loc").split(":")
            
        
        
        print("adding %s genome to diagram" % (record.name))
        seq_record = servers.bio_db.lookup(primary_id = record.name)
        
        """ Add the features from dictionary to the seq_record features """
        for location in locs:
            feature = SeqFeature(location = FeatureLocation(int(locs[location][0]), int(locs[location][1]), strand=1), type = location)
            seq_record.features.append(feature)
        """ Set up the Genome Diagram """
        max_len = max(max_len, len(seq_record))
        gd_track_for_features = gd_diagram.new_track(1, name = 
        seq_record.name, greytrack = True, start = 0, end = len(seq_record))
        gd_feature_set = gd_track_for_features.new_set()
        i = 0
            # Add Features
        print(seq_record.features)
        for feature in seq_record.features:
            print(feature)
            gd_feature_set.add_feature(feature, label = True, name
                               = str(i + 1), label_position = "start", label_size = 6, label_angle = 0)
            i+=1

    """ Draw and Write the Diagram to file """
    gd_diagram.draw(format="linear", pagesize = "A4", fragments = 0, start = 0, end = max_len)
    gd_diagram.write(output_path, "PNG")
    print("Genome Diagram has been added to file %s ", output_path)
    
def writeSeqToFile(ids):
    # Write Annotated Sequences to Genbank files to allow easy movement to Artemis
    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
    name = "GBTestSeq_"+phyloisland.randstring(5)
    out_path = "tmp/"+ name + ".gb"
    if len(query) > 1:
        # Multiple queries
        print("writing %s sequences to GenBank File" %(len(query)))
        seqlist = []
        for record in query.all():
            seq_record = servers.bio_db.lookup(primary_id = record.name)
            seqlist.append(seq_record)
        SeqIO.write(seqlist, out_path, "genbank")
    else:

        # Single Query
        seq_record = servers.bio_db.lookup(primary_id = query.name)
        SeqIO.write(seq_record, out_path, "genbank")
