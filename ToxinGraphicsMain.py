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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import glob 

def writeImageToFile(ids):

    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
    # Initiation of GenomeDiagram assets"
    name = "GenomeDiagram_"+phyloisland.randstring(5)
    gd_diagram = GenomeDiagram.Diagram(name)
    max_len = 0
    output_path = "tmp/"+name+".png"
    start = 0
    end = 0
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
            if feature.type in locs.keys():
                gd_feature_set.add_feature(feature, label = True, name
                                           = feature.type, color = region_colours[feature.type], label_position = "start", label_size = 6, label_angle = 0)
            elif feature.type == "CDS":
                 gd_feature_set.add_feature(feature)
                 
        """    For 'Zoomed' sections:     
        start = max(start, min(svals))
        if start > 1500:
            start -= 1000
        end = min(len(seq_record), max(endvals))
        if len(seq_record) - end > 1500:
            end += 1000
            """
            
    """ Draw and Write the Diagram to file """
    gd_diagram.draw(format="linear", pagesize = "A2", fragments = 0, start = start, end = end)
    gd_diagram.write(output_path, "PNG")
    print("Genome Diagram has been added to file %s ", output_path)
    
def writeSeqToFile(ids):
    # Write Annotated Sequences to Genbank files to allow easy movement to Artemis
    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
    print("writing sequences to GenBank File")
    for record in query.all():       
        """ Create a dictionary for key = feature type -> value = location """
        locs = {}
        colour_dict = {"a1":"255 165 0", "a2":"255 0 0", "a3":"255 255 0", "TcB":"0 0 255", "TcC":"255 0 255", "chitinase":"0 255 0"}
        # Prepare for literally the worst code in existence
        """ We have to pull the features from location data in database
        instead of directly from database because of current limitations """
        # Brute force if statements to pull from database
        # TODO - Could likely turn this into a for loop in some way
        if getattr(record, "a1_loc") != "":
            locs["a1"] = getattr(record, "a1_loc").split(":")
        if getattr(record, "a2 loc") != "":
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
        
    
        for location in locs:
            """ create and add features based on locations """
            color = {'color':colour_dict[location]}
            feature = SeqFeature(location = FeatureLocation(int(locs[location][0]), int(locs[location][1]), strand=1), type = "misc_feature", qualifiers = color)
            seq_record.features.append(feature)
            
        seq_record.seq = Seq(getattr(record, "sequence"), generic_dna)
        name = "GBTestSeq_"+record.name
        out_path = "tmp/"+ name + ".gb"
        SeqIO.write(seq_record, out_path, "genbank")



def writeHMMToImage(hmm_dict, reference, region, seq_record, species):
        # Currently taking region to be run simultaneously with hmmer, will need to manipulate some stuff for later
        # Initiation of GenomeDiagram assets"

        # Convert the length of the genome into the length of the translated genome.
        genome_length = round(len(seq_record) / 3)

        #Parameters that can be set

        # The highest number of tracks we'll attempt to add (if a feature overlaps on all tracks it'll just be added to the bottom track)
        # Set max_tracks to None to just add as many as we need to get no overlap
        max_tracks = None

        # The amount we want to allow a feature to overlap with a neighbouring feature
        overlap_amount = 0


        print("In hmm code")
        name = reference + "_" + species + "_GenomeDiagram"
        gd_diagram = GenomeDiagram.Diagram(name)
        max_len = 0
        output_path = reference + "/" + species + ".png"
        start = 0
        # For my work I was considering changing 'region1, 2, and 3' to a3, TcB, and TcC for convenience
        # Up to others though if I fully change that (is just a UI thing tbh)
        region_colours = {"a1": "orange", "a2": "red", "chitinase": "green", "a3": "yellow",
                          "TcB": "blue", "TcC": "magenta", "pore": "grey", "region1": "lightblue", "region2": "pink",
                          "region3": "purple", "region4": "black"}
        locs = {}
        strandd = 1
        for result in hmm_dict:
            i = 0
            for reg in result:
                """ Create a dictionary for key = feature type -> value = location """
                if "forward" in reg:
                    strandd = 1
                    locs[reg.split("/")[2] + phyloisland.randstring(5)] = result[reg].split(":")

                # If it's on the backwards strand we need to correct the postion - swap the start and end positions and set them relative to the 5' end
                # TODO: Actually this should be done in resultread, when we first read it in
                elif "backward" in reg:
                    strandd = -1
                    positions = result[reg].split(":")
                    print( "positions before split", positions)
                    positions[0], positions[1] = (len(seq_record) - int(positions[1])) , (len(seq_record) - int(positions[0]))
                    locs[reg.split("/")[2] + phyloisland.randstring(5)] = positions
                    print( "positions after split", positions)


                i += 1
                # Prepare for literally the worst code in existence
                """ We have to pull the features from location data in database
                instead of directly from database because of current limitations """
                # Brute force if statements to pull from database
                # TODO - Could likely turn this into a for loop in so

        """ Add the features from dictionary to the seq_record features """
        feature_locations = []
        # If want to add all features to image comment below line out
        seq_record.features = []
        i = 0
        for location in locs:
            """ Extract start and end values from each location, and add to independent lists """
            """ create and add features based on locations """
            feature = SeqFeature(
                location=FeatureLocation(int(locs[location][0]), int(locs[location][1]), strand=strandd),
                type=location[0:-5])
            seq_record.features.append(feature)

        """ Set up the Genome Diagram """
        max_len = max(max_len, len(seq_record))

        gd_track1 = gd_diagram.new_track(0, name=seq_record.name + " Track 1", greytrack=True, start=0,
                                         end=len(seq_record))
        gd_feature_set1 = gd_track1.new_set()

        for feature in seq_record.features:
            feature_locations.append((feature.location.start, feature.location.end))

        print("This many features")
        print(len(seq_record.features))

        total_tracks = 1
        # Dictionary to keep track of which locations are at which track
        forward_tracks = {1 : []}
        backward_tracks = {1 : []}

        for feature in seq_record.features:
            feature_added = False
            current_track = 1

            while current_track <= total_tracks and not feature_added:

                current_strand = feature.location.strand

                if current_strand == 1:
                    dict_track = forward_tracks
                elif current_strand == -1:
                    dict_track = backward_tracks


                overlap = False
                for loc in dict_track[current_track]:
                    if feature.location.start - overlap_amount in loc or feature.location.end + overlap_amount in loc:
                        overlap = True


                if overlap:

                    # If we've reached the highest track, make a new track
                    if current_track == total_tracks:
                        if max_tracks is not None and current_track == max_tracks:
                                # We've reached the highest track we want to try and add, so just add it to the bottom track (there will be overlap)
                                gd_feature_set1.add_feature(feature, label=True, name=feature.type, color=region_colours[feature.type], label_position='middle')
                                dict_track[1].append(feature.location)
                                feature_added = True
                        else:
                            # Add to this track
                            current_track += 1
                            total_tracks += 1
                            dict_track[current_track] = []
                            exec("gd_track" + str(current_track) + "= gd_diagram.new_track(0, name=seq_record.name + ' Track " + str(current_track) + "', greytrack=True, start=0, end=len(seq_record))")
                            exec("gd_feature_set" + str(current_track) + " = gd_track" + str(current_track) + ".new_set()")

                    else:
                        # Check on higher track
                        current_track += 1

                else:
                    # Add to this track
                    exec("gd_feature_set" + str(current_track) + ".add_feature(feature, label=True, name=feature.type, color=region_colours[feature.type], label_position='middle')")
                    if current_track in dict_track:
                        dict_track[current_track].append(feature.location)
                    else:
                        dict_track[current_track] = [feature.location]
                    feature_added = True


        for track in range(1, total_tracks):
            exec("gd_track"+ str(track) + ".add_set(gd_feature_set" + str(track) +")")
        gd_diagram.draw(format="linear", pagesize="A2", fragments=0, start=start, end=len(seq_record))
        gd_diagram.write(output_path, "PNG")
        print("Genome Diagram has been added to file " + output_path)


    
def writeHmmToSeq(hmm_dict, reference, region, seqrecord, species):
    name = species + "_sequence"
    output_path = reference +"/"+ name + ".gb"
    seqrecord.name = species[0:9]
    # Write Annotated Sequences to Genbank files to allow easy movement to Artemis
    print("writing sequences to GenBank File")       
    """ Create a dictionary for key = feature type -> value = location """
    locs = {}
    colour_dict = {"a1":"255 165 0", "a2":"255 0 0", "a3":"255 255 0", "TcB":"0 0 255", "TcC":"255 0 255",
                   "chitinase":"0 255 0", "region1":"0 255 255", "region2":"255 153 255", "region3":"204 0 102", "region4":"0 0 0"}
    for result in hmm_dict:
        i = 0
        for reg in result:       
            """ Create a dictionary for key = feature type -> value = location """
            if "forward" in reg:
                strandd = 1
            elif "backward" in reg:
                strandd = -1
            locs[reg.split("/")[2] + phyloisland.randstring(5)] = result[reg].split(":")
            i += 1
        
    print("adding %s genome to diagram" % (seqrecord.name))      
    
    seqrecord.features = []
    for location in locs:
        """ create and add features based on locations """
        color = {'color':colour_dict[location[0:-5]]}
        feature = SeqFeature(location = FeatureLocation(int(locs[location][0]), int(locs[location][1]), strand=strandd), type = "misc_feature", qualifiers = color)
        seqrecord.features.append(feature)
        
    sequence = str(seqrecord.seq)[2:-1]
    seqrecord.seq = Seq(sequence, generic_dna)
    SeqIO.write(seqrecord, output_path, "genbank")
