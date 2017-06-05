from flask import Flask, request, render_template, flash
from flask_sqlalchemy import SQLAlchemy
from flask_admin import Admin, form, BaseView, expose
from flask_admin.contrib.fileadmin import FileAdmin
from flask_admin.contrib.sqla import ModelView
from werkzeug.security import generate_password_hash, check_password_hash
from BioSQL import BioSeqDatabase
from forms import UploadForm
import os
import os.path as op
from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, SubmitField, TextAreaField, FileField
from wtforms.validators import DataRequired, Email, Length, ValidationError
import subMenu
from flask_uploads import UploadSet, configure_uploads, ALL
from Table import SortableTable, Item
from flask_admin.actions import action
import gettext
from Bio import Entrez, SeqIO, GenBank, AlignIO, pairwise2
from Bio.Align.Applications import MuscleCommandline
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
import requests
import subprocess
import sys
from Bio.SeqRecord import SeqRecord
from flask import Flask, make_response, request
from flask_admin.contrib.sqla import filters
from flask_admin.contrib.sqla.filters import BooleanEqualFilter
from flask_admin.contrib.sqla.filters import BaseSQLAFilter


# Setup temporary best guesses for what the A1 and A2 region should look like
yenA1 = "MDKYNNYSNVIKNKSSISPLLAAAAKIEPEITVLSSASKSNRSQYSQSLADTLLGLGYRSIFDIAKVSRQRFIKRHDESLLGNGAVIFDKAVSMANQVLQKYRKNRLEKSNSPLVPQTSSSTDASSESQTNKLPEYNQLFPEPWDNFCRPGAIEALDSPASYLLDLYKFIQSVELDGSNQARKLETRRADIPKLSLDNDALYKEVTALSIVNDVLSGSAREYIDQSGQADKAVNQILGDTHFPFTLPYSLPTQQINKGLGASNIELGTVIQRVDPQFSWNTTQEKYNQVLLAYTQLSSEQIALLSLPDVFTQNFLTQTELSAGYLSASTTEILAEKDLSRHGYIVKAADNIKGPTQLVEHSDASYDVIELTCTNQAKETITVKLRGENIITYQRTKARMVPFDNSSPFSRQLKLTFVAEDNPSLGNLDKGPYFANMDIYAAEWVRENVSSETMVSRPFLTMTYRIAIAKAGASLEELQPEADAFFINNFGLSAEDSSQLVKLVAFGDQTGSKAEEIESLLSCGENLPIVSPNVIFANPIFGSYFNDEPFPAPYHFGGVYINAHQRNAMTIIRAEGGREIQSLSNFRLERLNRFIRLQRWLDLPSHQLDLLLTSVMQADADNSQQEITEPVLKSLGLFRHLNLQYKITPEIFSSWLYQLTPFAVSGEIAFFDRIFNREQLFDQPFILDGGSFTYLDAKGSDAKSVKQLCAGLNISAVTFQFIAPLVQSALGLEAGTLVRSFEVVSSLYRLVSIPQTFGLSTEDGLILMNILTDEMGYLAKQPAFDDKQTQDKDFLSIILKMEALSAWLTKNNLTPASLALLLGVTRLAVVPTNNMVTFFKGIANGLSENVCLTTDDFQRQELEGADWWTLLSTNQVIDDMGLVLDIHPVWGKSDEEMLMEKIQSIGVSNDNNTLSIIVQILIQAKNAQENLLSQTISAEYGVERSVVPLQLRWLGSNVYSVLNQVLNNTPTDISSIVPKLSELTYSLLIYTQLINSLKLNKEFIFLRLTQPNWLGLTQPKLSTQLSLPEIYLITCYQDWVVNANKNEDSIHEYLEFANIKKTEAEKTLVDNSEKCAELLAEILAWDAGEILKAASLLGLNPPQATNVFEIDWIRRLQTLSEKTMISTEYLWQMGDLTENSEFSLKEGVGEAVMAALKAQGDSDNV"
yenA2 = "MSNSIEAKLQEDLRDALVDYYLGQIVPNSKDFINLRSTIKNVDDLYDHLLLDTQVSAKVITSRLSLVTQSVQQYINRIALNLEPGLSINQQEATDWEEFANRYGYWAANQQLRMFPEIYVDPTLRLTKTEFFFQLESALNQGKLTDDVAQKAVLGYLNNFEEVSNLEIIAGYQDGIDIENDKTYFVARTRMQPYRYFWRSLDASQRNSNSQELYPTAWSEWKAISVPLENVANGIVRPIMMDNRLYISWFEVAEEKDTDENGNIIVSGRYRTKIRLAHLGFDGIWSSGTTLREEVLAYQMEEMIAVVDRMEDEPRLALVAFKEMSENWDVVFSYICDSMLIESSNLPTTTHPPKPEDGDKGLSDLDDYGANLVWFYLHETANGGKAEYKQLILYPVIINRDWPIELDKTHQEGFGTVDDFTLNSNYTGDELSLYLQSSSTYKYDFSKSKNIIYGIWKEDANNNRCWLNYKLLTPEDYDPQINTTLVMCDKGDVNIITGFSLPNGGVDAGGKIKVTLRVGKKLRDKFQIKQFSQTQYLQFPEASSADVWYIGKQIRLNTLFAKELIGKASRSLDLVLSWETQNSRLEEAILGGAAELIDLDGANGIYFWELFFHMPFMVSWRFNVEQRYEDANRWIKYLFNPFECDDEPALLLGKPAYWNSRPLVDEPFTGYSLTQPSDPDAIAASDPIHYRKAVFNFLTKNIIDQGDMEYRKLQPSARTLARLSYSTASSLLGRRPDVQLTSFWQPLTLEDASYKTDSEIRAIEMQSQPLTFEPVVHDQTMSAVDNDIFMYPMNNELRGLWDRIENRIYNLRHNLTLDGKEINMDLYDSSISPRGLMKQRYQRVVTARNASKMNFKVPNYRFEPMLNRSKSGVETLIQFGSTLLSLLERKDSLSFDAYQMIQSGDLYRFSIDLQQQDIDINKASLEALQVSKQSAQDRYDHFKELYDENISSTEQKVIELQSQAANALLMAQGMRTAAAALDVIPNIYGLAVGGSHWGAPLNAAAEIIMIKYQADSSKSESLSVSESYRRRRQEWELQYKQAEWEVNSVEQQINLQNMQIKAANKRLEQVEAQQQQAMALLAYFSERFTNESLYTWLISQLSSLYLQAYDAVLSLCLSAEASLLYELNLGEQSFVGGGGWNDLYQGLMAGETLKLALMRMERVYVEQNSRRQEITKTISLKALLGESWPAELNKLKQKTPINFNLEEQIFVEDYQELYQRRIKSVSVSLPMLVGPYEDVCAQLTQTSSSYSTQADLKTVENMLTKRTFADTPHLVRSIQPNQQISLSTGVNDSGLFMLNFDDERFLPFEGSGVDSSWRLQFTNLKQNLDSLNDVILHVKYTAAIGSSTFSQGVRKILANINNDE"

# Create directory for file fields to use
file_path = op.join(op.dirname(__file__), 'files')
try:
    os.mkdir(file_path)
except OSError:
    pass

bio_server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="bioseqdb")
bio_db = bio_server["phylomain"]

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql://pi:@localhost/bioseqdb'
# app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql:///gabe'
app.config['SECRET_KEY'] = 'developmentkey'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = True

allfiles = UploadSet('all', ALL)
app.config['UPLOADS_ALL_DEST'] = 'static/uploads'
app.config['UPLOADED_ALL_DEST'] = 'static/uploads'
configure_uploads(app, allfiles)

db = SQLAlchemy(app)

# Create models
class File(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Unicode(64))
    path = db.Column(db.Unicode(128))

    def __unicode__(self):
        return self.name

class User(db.Model):
    __tablename__ = 'users'
    uid = db.Column(db.Integer, primary_key=True)
    firstname = db.Column(db.String(100))
    lastname = db.Column(db.String(100))
    studentno = db.Column(db.Integer)
    email = db.Column(db.String(120), unique=True)
    pwdhash = db.Column(db.String(54))

    def __init__(self, firstname="", lastname="", studentno="", email="", password=""):
        self.firstname = firstname.title()
        self.lastname = lastname.title()
        self.studentno = studentno
        self.email = email
        self.set_password(password)

    def set_password(self, password):
        self.pwdhash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.pwdhash, password)



# class BioEntry(db.Model):
#     __tablename__ = 'bioentry'
#     bioentry_id = db.Column(db.Integer, primary_key=True)
#     biodatabase_id = db.Column(db.Integer, unique = True)
#     taxon_id = db.Column(db.Integer)
#     name = db.Column(db.VARCHAR(40))
#     accession = db.Column(db.VARCHAR(128), unique=True)
#     identifier = db.Column(db.VARCHAR(40))
#     division = db.Column(db.VARCHAR(6))
#     description = db.Column(db.TEXT)
#     version = db.Column(db.SMALLINT, unique=True)
#
#     def __str__(self):
#         return str(self.bioentry_id)
#
# class BioEntryDBXRef(db.Model):
#     __tablename__ = 'bioentry_dbxref'
#     bioentry_id = db.Column(db.Integer, db.ForeignKey('bioentry.bioentry_id'), primary_key=True)
#     dbxref_id = db.Column(db.Integer)
#     rank = db.Column(db.Integer)
#     bioentry_name = db.Column(db.Column(db.VARCHAR(40)), db.ForeignKey('bioentry.name'))
#
#
#
#     relationship = db.relationship("BioEntry", backref="bioentry_dbxref")
#     # relationship = db.relationship("BioEntry", backref="bioentry_dbxref")
#
#
#     # form_ajax_refs = {
#     #     'relationship': {
#     #         'fields': ['name', 'description'],
#     #         'page_size': 10
#     #     }
#     # }


class SeqInstance(db.Model):
    __tablename__ = 'seqrecord'
    uid = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text)
    species = db.Column(db.String(255))
    strain = db.Column(db.String(255))
    description = db.Column(db.Text)
    a1 = db.Column(db.Text)
    a1_loc = db.Column(db.Text)
    a2 = db.Column(db.Text)
    a2_loc = db.Column(db.Text)

    overlap = db.Column(db.String(255))
    distance = db.Column(db.VARCHAR(255))
    sequence =  db.Column(db.Text)

    def __init__(self, name="", species="", strain="", description="", a1="", a1_loc="", a2="", a2_loc="", overlap="", distance="", sequence=""):
        self.name = name
        self.species = species
        self.strain = strain
        self.description = description
        self.a1 = a1
        self.a1_loc = a1_loc
        self.a2 = a2
        self.a2_loc = a2_loc
        self.overlap = overlap
        self.distance = distance
        self.sequence = sequence

class BioView(ModelView):
    column_sortable_list = ['name', 'division']




class FilterLastNameBrown(BaseSQLAFilter):
    def apply(self, query, value, alias=None):
        list = ['JGVH01000004.1', 'FMWJ01000001.1']

        return query.filter(self.column == "BX571867.1")


    def operation(self):
        return 'BX571867.1'

class FilterInAListMaybe(BaseSQLAFilter):

    def apply(self, query, value, alias=None):
        listup = ['JGVH01000004.1', 'FMWJ01000001.1']
        return query.filter(self.get_column(alias).in_(listup))

    def operation(self):
        return 'BX571867.1'


class FilterGetUnique(BaseSQLAFilter):
    list = ['JGVH01000004.1', 'FMWJ01000001.1']

    # def __init__(self, column, name, options=None, data_type=None):
    #     super(filters.FilterInList, self).__init__(column, name, options, data_type='select2-tags')
    #
    # def clean(self, value):
    #     return [v.strip() for v in value.split(',') if v.strip()]


    def apply(self, query, value, alias=None):
        return query.filter(self.get_column(query).in_(list))

    def operation(self):
        return 'in list'


class SeqInstanceView(ModelView):
    column_searchable_list = ['name', 'species', 'a1', 'a2', 'overlap']


    create_modal = True
    edit_modal = True

    def _a1description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        return model.a1[:15] + "..."

    def _a2description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        return model.a2[:15] + "..."

    def _seqdescription_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        return model.sequence[:15] + "..."
    column_formatters = {
        'a1': _a1description_formatter,
        'a2': _a2description_formatter,
        'sequence' : _seqdescription_formatter,
    }

    column_filters = ('name',
                      'species',
                      'strain',
                      'overlap',
                      FilterLastNameBrown(
                          SeqInstance.name, 'Da name is in'),
                      FilterInAListMaybe(
                          SeqInstance.name, 'Maybe this is in'),
                      # FilterLastNameBrown(column=SeqInstance.species, name="Get Unique Species Name"),
                      # filters.EqualFilter(column=SeqInstance.species, name="Photorhabdus luminescens"),
                      filters.FilterLike(SeqInstance.overlap, 'Overlaps', options=(('True', 'True'), ('False', 'False'))))
                      # BooleanEqualFilter(column=SeqInstance.species, name='Get unique species'))


    @action('getoverlap', 'Get Overlap')
    def action_getoverlap(self, ids):
        print (ids)
        try:
            query = SeqInstance.query.filter(SeqInstance.uid.in_(ids))

            count = 0
            for record in query.all():
                if record.a1 == "Not tested" or record.a2 == "Not tested":
                    pass
                else:
                    print (record.name, " is being checked")
                    distance = str(subMenu.getDistance(record.a1_loc, record.a2_loc))
                    if len(distance) > 1:
                        record.overlap = "False"
                        record.distance = distance
                    else:
                        record.overlap = "True"
                        record.distance = ""

                    db.session.add(record)
                    db.session.commit()
                    # print (subMenu.seqRecords)
                    # print (subMenu.seqRecords['CPYD01000004'].annotations.keys())

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')

    @action('checkA1', 'Check for A1 region')
    def action_checkA1(self, ids):
        print (ids)
        try:
            query = SeqInstance.query.filter(SeqInstance.uid.in_(ids))
            for record in query.all():
                print ("record name = ", record.name)
                seq_record = bio_db.lookup(primary_id=record.name)
                print ("seq record = ", seq_record.name)

                for feature in seq_record.features:
                    if 'translation' in feature.qualifiers:
                        print (feature.qualifiers['translation'][0])

                        alignment = pairwise2.align.globalms(yenA1, feature.qualifiers['translation'][0],  2, -1, -.5, -.1, score_only=True)
                        print (alignment)





            # query = SeqInstance.query.filter(SeqInstance.uid.in_(ids))

            # for record in query.all():
            #     if record.a1 != "Not tested":
            #         pass
            #     else:
            #         record = bio_db.lookup(id=)
            #         # print ('got here')
            #         # print (record)
            #         # print (record.name)
            #         # print (record.name.split(".")[0])
            #         # print (subMenu.seqRecords)
            #         # print (subMenu.seqRecords.get(record.name.split(".")[0]))
            #         # print (subMenu.seqRecords.get(record.name.split(".")[0]).features)
            #         #
            #         # print ("** KEYS **")
            #         # print (subMenu.seqRecords.get(record.name.split(".")[0]).annotations.keys())
            #         for feature in subMenu.seqDict.get(record.name.split(".")[0]).features:
            #             if 'gene' in feature.qualifiers and 'translation' in feature.qualifiers:
            #                 if "A1" in feature.qualifiers['gene'][0]:
            #                     print ("FOUND AN A1")
            #                     record.a1 = "Found"
            #                     db.session.add(record)
            #                     db.session.commit()

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')

    @action('checkA2', 'Check for A2 region')
    def action_checkA2(self, ids):
        print (ids)
        try:
            query = SeqInstance.query.filter(SeqInstance.uid.in_(ids))

            count = 0
            for record in query.all():
                if record.a1 != "Not tested":
                    pass
                else:
                    # print ('got here')
                    # print (record)
                    # print (record.name)
                    # print (record.name.split(".")[0])
                    # print (subMenu.seqRecords)
                    # print (subMenu.seqRecords.get(record.name.split(".")[0]))
                    # print (subMenu.seqRecords.get(record.name.split(".")[0]).features)
                    #
                    # print ("** KEYS **")
                    # print (subMenu.seqRecords.get(record.name.split(".")[0]).annotations.keys())
                    for feature in subMenu.seqDict.get(record.name.split(".")[0]).features:
                        if 'gene' in feature.qualifiers and 'translation' in feature.qualifiers:
                            if "A2" in feature.qualifiers['gene'][0]:
                                print ("FOUND AN A2")
                                record.a1 = "Found"
                                db.session.add(record)
                                db.session.commit()

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')

    @action('buildprofileA2', 'Build profile based on A2')
    def action_build_profile_a2(self, ids):
        print (ids)
        try:
            query = SeqInstance.query.filter(SeqInstance.uid.in_(ids))
            align_list = []
            for record in query.all():
                if record.a2 == "Not tested":
                    pass
                else:
                    print (record.a2)
                    print (record.name)
                    align_record = SeqRecord(Seq(record.a2, generic_protein), id=str(record.name) + "_" + "A2")
                    # newSeq = Seq(record.a2)
                    # align_record = Seq(newSeq, )
                    # align_record = SeqRecord(Seq(record.a2, generic_protein), id="gabe")
                    align_list.append(align_record)

                    # Write sequences to FASTA

            SeqIO.write(align_list, "align.fasta", "fasta")
            muscle_cline = MuscleCommandline(input="align.fasta")
            child = subprocess.Popen(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                         universal_newlines=True, shell=(sys.platform != "win32"))
            alignment = AlignIO.read(child.stdout, "fasta")
            AlignIO.write(alignment, "align.aln", "fasta")
            file = subprocess.call(["hmmbuild", "profile.hmm", "align.aln"])

            # file_contents = file.stream.read().decode("utf-8")




        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')

    @action('showyen', 'Show Yen')
    def action_build_profile_a2(self, ids):
        print (ids)
        try:
            print ('trying')
            query = self.session.query(self.model).filter(self.model.species == "Yersinia nurmii")
            render_template("index.html")
            return query

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')

class FilterGetBrown(BaseSQLAFilter):
    def apply(self, query, value, alias=None):
        list = ['JGVH01000004.1', 'FMWJ01000001.1' ]
        if value == '1':
            return query.filter(self.column=="Brown")
        else:
            return query.filter(self.column != "Brown")

        def operation(self):
            return 'is Brown'





class UploadForm(FlaskForm):
    name = StringField('What ID should we give this feature?', validators=[DataRequired("Not completed")])
    file = FileField('Upload the FASTA file', validators=[DataRequired("Not completed")] )

    upload_submit = SubmitField("Upload files")

class SelectedView(ModelView):
    def get_query(self):
        return self.session.query(self.model).filter(self.model.species=="Yersinia nurmii")

    # def get_count_query(self):
    #     return self.session.query('*').filter(self.model.species=="Yersinia nurmii")

class FileUpload(FileAdmin):
    can_mkdir = False
    can_delete = False
    can_delete_dirs = False

class UserView(ModelView):
    can_create = False
    column_exclude_list = ["pwdhash"]
    @action('approve', 'Approve', 'Are you sure you want to approve selected users?')
    def action_approve(self, ids):
        print (ids)
        try:
            query = User.query.filter(User.uid.in_(ids))

            count = 0
            for user in query.all():
                print (user)
                if user.firstname == 'Gabe':
                    count += 1

                print ("Count was ", count)
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')


class MyForm(ModelView):
    def scaffold_form(self):
        form_class = super(UserView, self).scaffold_form
        form_class.extra = StringField('Extra')
        return form_class


class UploadView(BaseView):


    # @app.route("/upload", methods = ['GET', 'POST'])
    @expose("/", methods =('GET', 'POST'))
    def upload(self):
        form = UploadForm()
        sort = request.args.get('sort', 'id')
        reverse = (request.args.get('direction', 'asc') == 'desc')

        if request.method == 'POST':
            # region = request.form['name']

            # Get the file and ID to name these sequences
            region = form.name.data
            filename = allfiles.save(request.files['file'])
            print (filename, region, "##")

            # Create the initial seqRecords
            subMenu.seqDict = {}
            subMenu.unmappable = []
            subMenu.defaultValue("static/uploads/"  + filename, region)
            print ("@@@@@@@@@@@@@@Couldn't map@@@@@@@@@@@@@" , len(subMenu.unmappable))

            records = subMenu.seqDict




            for record in records:
                current = records[record]
                # for feat in current.features:
                #     print ("**", feat)
                print (current)
                name = current.id
                species = current.annotations['species']
                strain = current.annotations['source']
                description = current.description
                a1 = current.annotations["A1"] if "A1" in current.annotations.keys() else "Not tested"
                a1_loc = current.annotations["A1_location"] if "A1_location" in current.annotations.keys() else "Not tested"
                a2 = current.annotations["A2"] if "A2" in current.annotations.keys() else "Not tested"
                a2_loc = current.annotations["A2_location"] if "A2_location" in current.annotations.keys() else "Not tested"
                overlap = current.annotations["Overlap"] if "Overlap" in current.annotations.keys() else "Not tested"
                distance = "Not tested"
                sequence = str(current.seq)

                entry = SeqInstance(name, species, strain, description, a1, a1_loc, a2, a2_loc, overlap, distance, sequence)

                check = SeqInstance.query.filter_by(name=name).first()

                # check = db.session.query(db.exists().where(SeqInstance.name== name)).scalar()

                print ("CHECKING BABY")
                print ("Name is" , name)
                # print (current)
                print (check)

                if (check):
                    print('Name we are checking on is ', name)

                    print ("A1  ", check.a1)
                    print ("A2  ", check.a2)
                    print ("region is ", region)
                    if region in current.annotations.keys():
                        print ("YEAH BABY")
                    else:
                        print ("NAH BOY")
                    setattr(check, region.lower(), current.annotations[region] if region in current.annotations.keys() else "Not tested")
                    setattr(check, region.lower() + "_loc", current.annotations[region + "_location"] if region + "_location" in current.annotations.keys() else "Not tested")
                    # print ('got here baby')

                    # check = subMenu.addToRecord(check, current, region)
                    # print ("A1  ", check.a1)
                    # print ("A2  ", check.a2)

                    db.session.add(check)
                else:
                    print('Name apparently not there is ', name)


                    items = subMenu.seqDict.items()
                    print ("Adding to db")
                    seqList = []

                    db.session.add(entry)
                    # for k, v in items:
                    #     print (record)
                    #     print (v)
                    #     seqList.append(v)

                    seqList.append(current)
                    print ("Adding to bio_db")
                    count = bio_db.load(seqList)
                    print ("Making the bio commit")
                    bio_server.commit()

                    # print('Loaded %d sequences' % count)

                ("Making the seqrecord commit")
                db.session.commit()

            return self.render("upload_admin.html", form=form, records = records)
        return self.render("upload_admin.html", form=form)



admin = Admin(app, name="Phylo Island", template_mode="bootstrap3")
admin.add_view(UploadView(name='Upload', endpoint='upload_admin'))

admin.add_view(SeqInstanceView(SeqInstance, db.session, endpoint="seq_view")) # working version



# admin.add_view(BioView(BioEntry, db.session))
# admin.add_view(UserView(User, db.session))
# admin.add_view(BioEntryDBXRefView(BioEntryDBXRef, db.session))


# admin.add_view(SeqRecordView(SeqRecord, db.session))



# path = op.join(op.dirname(__file__), 'static') # path for static view
# admin.add_view(FileAdmin(path, '/static/', name='Static Files')) # static view


# admin.add_view(SelectedView(SeqInstance, db.session, endpoint="selected_view")) # selected view


# admin.add_view(FileView(File, db.session))
# admin.add_view(FileView(File, db.session))
# admin.add_view(UserView(User, db.session, name='User'))
# admin.add_view(AnalyticsView(name='Analytics', endpoint='analytics'))


@app.route('/')
def index():
    form = UploadForm()
    # return render_template("index.html", form=form)
    return '<a href="/admin/">Click me to get to Admin!</a>'

# @app.route("/upload", methods = ['GET', 'POST'])
# def upload():
#
#     form = UploadForm()
#     return render_template("admin/upload_admin.html", form=form)


@app.route("/admin", methods = ['GET', 'POST'])
def admin():
    form = UploadForm()
    return render_template("admin/index.html", form=form)

    if request.method == 'POST':
        record = bio_db.lookup(accession="CPYD01000004")
        print(record)
    else:
        return render_template("admin/index.html", form=form)





if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0')

