from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from reportlab.lib import colors

# Read in the test GenBank file
record = SeqIO.read("./testing/test.gb", "genbank")

# Create a new GenomeDiagram
gd_diagram = GenomeDiagram.Diagram(record.id)

# Create two tracks to store the features
gd_track1 = gd_diagram.new_track(1, name="Track 1")
gd_track2 = gd_diagram.new_track(2, name="Track 2")

# Create the feature sets
gd_feature_set1 = gd_track1.new_set()
gd_feature_set2 = gd_track2.new_set()

feature_list = []
# Create all the feature sets
first_feature = SeqFeature(FeatureLocation(10000, 50000, strand=+1))
feature_list.append(first_feature)
second_feature = SeqFeature(FeatureLocation(60000, 80000, strand=+1))
feature_list.append(second_feature)
third_feature = SeqFeature(FeatureLocation(90000, 120000, strand=+1))
feature_list.append(third_feature)
fourth_feature = SeqFeature(FeatureLocation(110000, 150000, strand=+1))
feature_list.append(fourth_feature)
overlapping_feature = SeqFeature(FeatureLocation(10000, 70000, strand=+1))
feature_list.append(overlapping_feature)

gd_feature_set_copy = gd_feature_set1

i = 0
feature_locations = []
for feature in feature_list:
    feature_locations.append((feature.location.start, feature.location.end))

for feature in feature_list:
    if i == 0:
        gd_feature_set1.add_feature(feature_list[i], color = colors.green)
    while i +1 < len(feature_locations):
        overlap = feature.location.start in feature_locations[i + 1] or feature.location.end in feature_locations[i+ 1]
        i +=1
        print(overlap)
        if overlap:
            gd_feature_set2.add_feature(feature_list[i], color = colors.red)
        else:
            gd_feature_set1.add_feature(feature_list[i], color = colors.green, label = True)



    
# Now we're adding all the features to the first feature set
#gd_feature_set1.add_feature(first_feature, name="First feature", color=colors.green, label_color=colors.green, label_position='middle', label=True)
#gd_feature_set1.add_feature(second_feature, name="Second feature", color=colors.blue, label_color=colors.blue, label_position='middle', label=True)
#gd_feature_set1.add_feature(third_feature, name="Third feature", color=colors.orange, label_color=colors.orange, label_position='middle', label=True)
#gd_feature_set1.add_feature(fourth_feature, name="Fourth feature", color=colors.brown, label_color=colors.brown, label_position='middle', label=True)

# Note that in the actual code we'd want to attempt to add this feature and see that it would overlap - and so create
# a new track and add the feature to this track instead
#gd_feature_set2.add_feature(overlapping_feature, name="Overlapping feature", color=colors.red,  label_color=colors.red, label_position='middle', label=True)


# Add the feature sets to the track
gd_track1.add_set(gd_feature_set1)
gd_track2.add_set(gd_feature_set2)


# Draw the diagram
gd_diagram.draw(format="linear", pagesize='A4', fragments=4,
                start=0, end=len(record))

# Write it to file
gd_diagram.write("./testing/CP010029.pdf", "PDF")
