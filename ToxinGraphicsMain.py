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
    #start = 0
    #end = 0
    # For my work I was considering changing 'region1, 2, and 3' to a3, TcB, and TcC for convenience
    # Up to others though if I fully change that (is just a UI thing tbh)
    region_colours = {"a1":"orange", "a2":"red", "chi":"green", "a3":"yellow",
                      "TcB":"blue", "TcC":"magenta", "pore":"grey", "region4":"pink"}
    for record in query.all():
        
        """ Create a dictionary for key = feature type -> value = location """
        locs = {}
        # Prepare for literally the worst code in existence
        """ We have to pull the features from location data in database
        instead of directly from database because of current limitations """
        # Brute force if statements to pull from database
        # TODO - Could likely turn this into a for loop in some way
        if getattr(record, "a1_loc") != "":
            locs["a1"] = getattr(record, "a1_loc").split(":")
        if getattr(record, "a2_loc") != "":
            locs["a2"] = getattr(record, "a2_loc").split(":")
        if getattr(record, "pore_loc") != "":
            locs["pore"] = getattr(record, "pore_loc").split(":")
        if getattr(record, "chitinase_loc") != "":
            locs["chi"] = getattr(record, "chitinase_loc").split(":")
        if getattr(record, "region1_loc") != "":
            locs["a3"] = getattr(record, "region1_loc").split(":")
        if getattr(record, "region2_loc") != "":
            locs["TcB"] = getattr(record, "region2_loc").split(":")
        if getattr(record, "region3_loc") != "":
            locs["TcC"] = getattr(record, "region3_loc").split(":")
        if getattr(record, "region1_loc") != "":
            locs["region4"] = getattr(record, "region4_loc").split(":")
        
        
        print("adding %s genome to diagram" % (record.name))
        seq_record = servers.bio_db.lookup(primary_id = record.name)
        
        """ Add the features from dictionary to the seq_record features """
        svals = []
        endvals = []        
        for location in locs:
            """ Extract start and end values ffrom each location, and add to independent lists """
            sval = int(locs[location][0])
            endval = int(locs[location][1])
            svals.append(sval)
            endvals.append(endval)
            """ create and add features based on locations """
            feature = SeqFeature(location = FeatureLocation(int(locs[location][0]), int(locs[location][1]), strand=1), type = location)
            seq_record.features.append(feature)

        """ Set up the Genome Diagram """
        max_len = max(max_len, len(seq_record))

        gd_track_for_features = gd_diagram.new_track(1, name = 
        seq_record.name, greytrack = True, start = 0, end = len(seq_record))
        gd_feature_set = gd_track_for_features.new_set()
            # Add Features
        print(seq_record.features)
        for feature in seq_record.features:
            print(feature)
            gd_feature_set.add_feature(feature, label = True, name
                               = feature.type, color = region_colours[feature.type], label_position = "start", label_size = 6, label_angle = 0)
        """start = max(start, min(svals))
        if start > 1500:
            start -= 1000
        end = min(len(seq_record), max(endvals))
        if len(seq_record) - end > 1500:
            end += 1000
            """
    """ Draw and Write the Diagram to file """
    gd_diagram.draw(format="linear", pagesize = "A4", fragments = 1, start = 0, end = max_len)
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
