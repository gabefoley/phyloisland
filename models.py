from servers import *
from werkzeug.datastructures import FileStorage
from wtforms import ValidationError, fields
from wtforms.validators import required
from wtforms.widgets import HTMLString, html_params, FileInput
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, FileField, SelectField, BooleanField, validators
from wtforms.validators import DataRequired
from gettext import gettext

# Create models

class SequenceRecords(db.Model):
    __tablename__ = 'seqrecord'
    uid = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text)
    species = db.Column(db.String(255))
    description = db.Column(db.Text)
    sequence = db.Column(db.Text)
    a1_ref = db.Column(db.Boolean)
    a2_ref = db.Column(db.Boolean)

    def __init__(self, name="", species="", description="", sequence="", a1_ref=0, a2_ref=0):
        self.name = name
        self.species = species
        self.description = description
        self.sequence = sequence
        self.a1_ref = a1_ref
        self.a2_ref = a2_ref


class GenomeRecords(db.Model):
    __tablename__ = 'genomerecord'
    uid = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text)
    species = db.Column(db.String(255))
    strain = db.Column(db.String(255))
    description = db.Column(db.Text)
    a1 = db.Column(db.Text)
    a1_length = db.Column(db.Integer)
    a1_loc = db.Column(db.Text)
    a2 = db.Column(db.Text)
    a2_length = db.Column(db.Integer)
    a2_loc = db.Column(db.Text)

    overlap = db.Column(db.String(255))
    distance = db.Column(db.VARCHAR(255))
    sequence = db.Column(db.Text)
    a1_ref = db.Column(db.Boolean)
    a2_ref = db.Column(db.Boolean)


    def __init__(self, name="", species="", strain="", description="", a1="", a1_length="", a1_loc="", a2="", a2_length="", a2_loc="", overlap="",
                 distance="", sequence="", a1_ref=0, a2_ref=0):
        self.name = name
        self.species = species
        self.strain = strain
        self.description = description
        self.a1 = a1
        self.a1_length = a1_length
        self.a1_loc = a1_loc
        self.a2 = a2
        self.a2_length = a2_length
        self.a2_loc = a2_loc
        self.overlap = overlap
        self.distance = distance
        self.sequence = sequence
        self.a1_ref = a1_ref
        self.a2_ref = a2_ref


class BlobMixin(object):
    mimetype = db.Column(db.Unicode(length=255), nullable=False)
    filename = db.Column(db.Unicode(length=255), nullable=False)
    profile = db.Column(db.BLOB, nullable=False)
    size = db.Column(db.Integer, nullable=False)

    def __init__(self, mimetype, filename, profile, size):
        self.mimetype = mimetype
        self.filename = filename
        self.profile = profile
        self.size = size


class Profile(db.Model, BlobMixin):
    __tablename__ = 'profile_blob'

    uid = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Unicode(length=255), nullable=False, unique=True)
    a1_profile_ref = db.Column(db.Boolean)
    a2_profile_ref = db.Column(db.Boolean)

    def __init__(self, name = "", blobMix = "", a1_ref=0, a2_ref=0):
        self.name = name
        self.a1_ref = a1_ref
        self.a2_ref = a2_ref
        self.blobMix = blobMix

    def set_blobMix(self, blobMix):
        self.blobMix = blobMix
        self.mimetype = blobMix.mimetype
        self.filename = blobMix.filename
        self.profile = blobMix.profile
        self.size = blobMix.size



    def __unicode__(self):
        return u"name : {name}; filename : {filename})".format(name=self.name, filename=self.filename)


class BlobUploadField(fields.StringField):

    widget = FileInput()

    def __init__(self, label=None, allowed_extensions=None, size_field=None, filename_field=None, mimetype_field=None, **kwargs):

        self.allowed_extensions = allowed_extensions
        self.size_field = size_field
        self.filename_field = filename_field
        self.mimetype_field = mimetype_field
        validators = [required()]

        super(BlobUploadField, self).__init__(label, validators, **kwargs)

    def is_file_allowed(self, filename):
        """
            Check if file extension is allowed.

            :param filename:
                File name to check
        """
        if not self.allowed_extensions:
            return True

        return ('.' in filename and
                filename.rsplit('.', 1)[1].lower() in
                map(lambda x: x.lower(), self.allowed_extensions))

    def _is_uploaded_file(self, data):
        return (data and isinstance(data, FileStorage) and data.filename)

    def pre_validate(self, form):
        super(BlobUploadField, self).pre_validate(form)
        if self._is_uploaded_file(self.data) and not self.is_file_allowed(self.data.filename):
            raise ValidationError(gettext('Invalid file extension'))

    def process_formdata(self, valuelist):
        if valuelist:
            data = valuelist[0]
            self.data = data

    def populate_obj(self, obj, name):


        if self._is_uploaded_file(self.data):

            _profile = self.data.read()



            setattr(obj, name, _profile)


            # setattr(obj, self.profile, _blob)


            if self.size_field:
                setattr(obj, self.size_field, len(_profile))

            if self.filename_field:
                setattr(obj, self.filename_field, self.data.filename)

            if self.mimetype_field:
                setattr(obj, self.mimetype_field, self.data.content_type)

# Form for uploading files
class UploadForm(FlaskForm):
    file = FileField('Upload the file that contains the information we will map to the genome records.', [validators.DataRequired()])
    type = SelectField('What type of file is this?', [validators.DataRequired()], choices = [("protein", "FASTA (amino acids)"), ("nucleotide", "FASTA (nucleotides)"), ("species", "Species list"), ("genome", "Genome ID list")])
    add_sequence = BooleanField("Add sequences to sequence database?", default="checked")


    upload_submit = SubmitField("Upload file")


