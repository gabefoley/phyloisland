
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO, pairwise2
import servers
import models
from flask import flash
import re

def checkFeature(ids, reference, recordName, recordLocation):
    best_align = 0
    best_location = None
    best_seq = None
    query = models.SequenceRecords.query.filter(models.SequenceRecords.uid.in_(ids))
    for record in query.all():
        print(record.name)
        seq_record = servers.bio_db.lookup(primary_id=record.name)

        print(seq_record.description)
        print(seq_record)

        for feature in seq_record.features:
            print(feature)
            if 'translation' in feature.qualifiers:

                alignment = pairwise2.align.globalms(reference, feature.qualifiers['translation'][0],
                                                     2, -1, -1, -.5, score_only=True)
                if alignment > best_align:
                    best_align = alignment
                    best_seq = feature.qualifiers['translation'][0]
                    best_location = feature.location

        if (best_seq):
            print("best align is ", best_align)
            flash("Found an A1 region in %s" % record.name)
            location = re.search(r"\d*:\d*", str(best_location))

            setattr(record, recordName, best_seq)
            print('location is')
            print(location)
            if location:
                print(location.group(0))
                setattr(record, recordLocation, location.group(0))
                servers.db.session.add(record)
                servers.db.session.commit()

        # alignment = pairwise2.align.globalms(yenA1, best_seq,
        #                                      2, -1, -1, -.5)