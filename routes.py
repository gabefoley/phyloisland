
from typing import Any
from os.path import join
from flask import flash
from flask_sqlalchemy import SQLAlchemy
from flask_admin import Admin, BaseView, expose
from flask_admin.contrib.sqla import ModelView
from werkzeug.security import generate_password_hash, check_password_hash
from BioSQL import BioSeqDatabase
import os
import os.path as op
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, FileField
from wtforms.validators import DataRequired
import phyloisland
from flask_uploads import UploadSet, configure_uploads, ALL
from flask_admin.actions import action
import gettext
from Bio import SeqIO, AlignIO, pairwise2
from Bio.Align.Applications import MuscleCommandline
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
import subprocess
import sys
from Bio.SeqRecord import SeqRecord
from flask import Flask, request
from flask_admin.contrib.sqla import filters
from flask_admin.contrib.sqla.filters import BaseSQLAFilter
from flask_admin import AdminIndexView
import re



# Setup temporary best guesses for what the A1 and A2 region should look like
yenA1 = "MDKYNNYSNVIKNKSSISPLLAAAAKIEPEITVLSSASKSNRSQYSQSLADTLLGLGYRSIFDIAKVSRQRFIKRHDESLLGNGAVIFDKAVSMANQVLQKYRKNRL" \
        "EKSNSPLVPQTSSSTDASSESQTNKLPEYNQLFPEPWDNFCRPGAIEALDSPASYLLDLYKFIQSVELDGSNQARKLETRRADIPKLSLDNDALYKEVTALSIVNDV" \
        "LSGSAREYIDQSGQADKAVNQILGDTHFPFTLPYSLPTQQINKGLGASNIELGTVIQRVDPQFSWNTTQEKYNQVLLAYTQLSSEQIALLSLPDVFTQNFLTQTELS" \
        "AGYLSASTTEILAEKDLSRHGYIVKAADNIKGPTQLVEHSDASYDVIELTCTNQAKETITVKLRGENIITYQRTKARMVPFDNSSPFSRQLKLTFVAEDNPSLGNLD" \
        "KGPYFANMDIYAAEWVRENVSSETMVSRPFLTMTYRIAIAKAGASLEELQPEADAFFINNFGLSAEDSSQLVKLVAFGDQTGSKAEEIESLLSCGENLPIVSPNVIF" \
        "ANPIFGSYFNDEPFPAPYHFGGVYINAHQRNAMTIIRAEGGREIQSLSNFRLERLNRFIRLQRWLDLPSHQLDLLLTSVMQADADNSQQEITEPVLKSLGLFRHLNL" \
        "QYKITPEIFSSWLYQLTPFAVSGEIAFFDRIFNREQLFDQPFILDGGSFTYLDAKGSDAKSVKQLCAGLNISAVTFQFIAPLVQSALGLEAGTLVRSFEVVSSLYRL" \
        "VSIPQTFGLSTEDGLILMNILTDEMGYLAKQPAFDDKQTQDKDFLSIILKMEALSAWLTKNNLTPASLALLLGVTRLAVVPTNNMVTFFKGIANGLSENVCLTTDDF" \
        "QRQELEGADWWTLLSTNQVIDDMGLVLDIHPVWGKSDEEMLMEKIQSIGVSNDNNTLSIIVQILIQAKNAQENLLSQTISAEYGVERSVVPLQLRWLGSNVYSVLNQ" \
        "VLNNTPTDISSIVPKLSELTYSLLIYTQLINSLKLNKEFIFLRLTQPNWLGLTQPKLSTQLSLPEIYLITCYQDWVVNANKNEDSIHEYLEFANIKKTEAEKTLVDN" \
        "SEKCAELLAEILAWDAGEILKAASLLGLNPPQATNVFEIDWIRRLQTLSEKTMISTEYLWQMGDLTENSEFSLKEGVGEAVMAALKAQGDSDNV"

yenA2 = "MSNSIEAKLQEDLRDALVDYYLGQIVPNSKDFINLRSTIKNVDDLYDHLLLDTQVSAKVITSRLSLVTQSVQQYINRIALNLEPGLSINQQEATDWEEFANRYGYWA" \
        "ANQQLRMFPEIYVDPTLRLTKTEFFFQLESALNQGKLTDDVAQKAVLGYLNNFEEVSNLEIIAGYQDGIDIENDKTYFVARTRMQPYRYFWRSLDASQRNSNSQELY" \
        "PTAWSEWKAISVPLENVANGIVRPIMMDNRLYISWFEVAEEKDTDENGNIIVSGRYRTKIRLAHLGFDGIWSSGTTLREEVLAYQMEEMIAVVDRMEDEPRLALVAF" \
        "KEMSENWDVVFSYICDSMLIESSNLPTTTHPPKPEDGDKGLSDLDDYGANLVWFYLHETANGGKAEYKQLILYPVIINRDWPIELDKTHQEGFGTVDDFTLNSNYTG" \
        "DELSLYLQSSSTYKYDFSKSKNIIYGIWKEDANNNRCWLNYKLLTPEDYDPQINTTLVMCDKGDVNIITGFSLPNGGVDAGGKIKVTLRVGKKLRDKFQIKQFSQTQ" \
        "YLQFPEASSADVWYIGKQIRLNTLFAKELIGKASRSLDLVLSWETQNSRLEEAILGGAAELIDLDGANGIYFWELFFHMPFMVSWRFNVEQRYEDANRWIKYLFNPF" \
        "ECDDEPALLLGKPAYWNSRPLVDEPFTGYSLTQPSDPDAIAASDPIHYRKAVFNFLTKNIIDQGDMEYRKLQPSARTLARLSYSTASSLLGRRPDVQLTSFWQPLTL" \
        "EDASYKTDSEIRAIEMQSQPLTFEPVVHDQTMSAVDNDIFMYPMNNELRGLWDRIENRIYNLRHNLTLDGKEINMDLYDSSISPRGLMKQRYQRVVTARNASKMNFK" \
        "VPNYRFEPMLNRSKSGVETLIQFGSTLLSLLERKDSLSFDAYQMIQSGDLYRFSIDLQQQDIDINKASLEALQVSKQSAQDRYDHFKELYDENISSTEQKVIELQSQ" \
        "AANALLMAQGMRTAAAALDVIPNIYGLAVGGSHWGAPLNAAAEIIMIKYQADSSKSESLSVSESYRRRRQEWELQYKQAEWEVNSVEQQINLQNMQIKAANKRLEQV" \
        "EAQQQQAMALLAYFSERFTNESLYTWLISQLSSLYLQAYDAVLSLCLSAEASLLYELNLGEQSFVGGGGWNDLYQGLMAGETLKLALMRMERVYVEQNSRRQEITKT" \
        "ISLKALLGESWPAELNKLKQKTPINFNLEEQIFVEDYQELYQRRIKSVSVSLPMLVGPYEDVCAQLTQTSSSYSTQADLKTVENMLTKRTFADTPHLVRSIQPNQQI" \
        "SLSTGVNDSGLFMLNFDDERFLPFEGSGVDSSWRLQFTNLKQNLDSLNDVILHVKYTAAIGSSTFSQGVRKILANINNDE"

# Create directory for file fields to use
file_path = op.join(op.dirname(__file__), 'files')
try:
    os.mkdir(file_path)
except OSError:
    pass
BASE_ROUTE = '/phyloisland'

def local(route: str) -> str:
    if BASE_ROUTE == '/':
        return route
    else:
        return join(BASE_ROUTE, route[1:])

def local_url_for(*args, **kwargs) -> str:
    new_url = local(url_for(*args, **kwargs))
    if new_url.count(BASE_ROUTE[1:]) == 1:
        fixed_url = new_url
        return new_url
    else:
        fixed_url = '/'.join(new_url.split('/')[2:])
    assert fixed_url.count(BASE_ROUTE[1:]) == 1, fixed_url
    return fixed_url

def local_redirect(*args, **kwargs) -> Any:
    return redirect(local_url_for(*args, **kwargs))

bio_server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="fishtank")
# bio_db = bio_server["fishface"]
bio_db = bio_server["fish"]


application = Flask(__name__)
application.config['SQLALCHEMY_DATABASE_URI'] = 'mysql://pi:@localhost/fishtank'
application.config['SECRET_KEY'] = 'developmentkey'
application.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = True

allfiles = UploadSet('all', ALL)
application.config['UPLOADS_ALL_DEST'] = 'static/uploads'
application.config['UPLOADED_ALL_DEST'] = 'static/uploads'
configure_uploads(application, allfiles)

db = SQLAlchemy(application)


# Create models
class SequenceRecords(db.Model):
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
    sequence = db.Column(db.Text)

    def __init__(self, name="", species="", strain="", description="", a1="", a1_loc="", a2="", a2_loc="", overlap="",
                 distance="", sequence=""):
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


class Profiles(db.Model):
    __tablename__ = 'profile'
    uid = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text)
    profile = db.Column(db.BLOB)

    def __init__(self, name="", profile=""):
        self.name = name
        self.profile = profile


class FilterInAListMaybe(BaseSQLAFilter):

    def apply(self, query, value, alias=None):
        listup = ['JGVH01000004.1', 'FMWJ01000001.1', 'CPYD01000004.1', 'DQ400808.1']
        return query.filter(self.get_column(alias).in_(listup))

    def operation(self):
        return 'Yes'




class FilterGetUnique(BaseSQLAFilter):
    list = ['JGVH01000004.1', 'FMWJ01000001.1', 'CPYD01000004.1', 'DQ400808.1' ]

    # def __init__(self, column, name, options=None, data_type=None):
    #     super(filters.FilterInList, self).__init__(column, name, options, data_type='select2-tags')
    #
    # def clean(self, value):
    #     return [v.strip() for v in value.split(',') if v.strip()]

    def apply(self, query, value, alias=None):
        return query.filter(self.get_column(query).in_(list))

    def operation(self):
        return 'in list'


class SequenceRecordsView(ModelView):
    column_searchable_list = ['name', 'species', 'a1', 'a2', 'overlap']
    create_modal = True
    edit_modal = True
    can_create = False
    can_view_details = True

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
        'sequence': _seqdescription_formatter,
    }

    column_filters = ('name',
                      'species',
                      'strain',
                      'description',
                      'a1',
                      'a2',
                      filters.FilterLike(SequenceRecords.overlap, 'Overlap',
                                         options=(('True', 'True'), ('False', 'False'))),
                      'sequence',
                      FilterInAListMaybe(
                          SequenceRecords.name, 'Get unique species'))

    # Delete the record from the BioSQL database as well
    def after_model_delete(self, model):
        print ("model gone")
        print (model.uid)
        print (model.name)
        del bio_db[model.uid]
        bio_server.commit()





    @action('do_something', 'Do Something')
    def do_something(self, ids):
        for id in ids:
            del bio_db[str(id)]
        bio_server.commit()



    @action('c_getoverlap', 'Get Overlap')
    def action_getoverlap(self, ids):
        try:
            query = SequenceRecords.query.filter(SequenceRecords.uid.in_(ids))
            for record in query.all():
                if record.a1 == "Not tested" or record.a2 == "Not tested":
                    pass
                else:
                    distance = str(phyloisland.getDistance(record.a1_loc, record.a2_loc))
                    if len(distance) > 1:
                        record.overlap = "False"
                        record.distance = distance
                    else:
                        record.overlap = "True"
                        record.distance = ""

                    db.session.add(record)
                    db.session.commit()

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')

    @action('a_check_a1', 'Check for A1 region')
    def action_check_a1(self, ids):
        try:
            best_align = 0
            query = SequenceRecords.query.filter(SequenceRecords.uid.in_(ids))
            for record in query.all():
                print(record.name)
                seq_record = bio_db.lookup(primary_id=record.name)

                for feature in seq_record.features:
                    if 'translation' in feature.qualifiers:
                        alignment = pairwise2.align.globalms(yenA1, feature.qualifiers['translation'][0],
                                                             2, -1, -.5, -.1, score_only=True)
                        if alignment > best_align:
                            best_align = alignment
                            best_seq = feature.qualifiers['translation'][0]
                            best_location = feature.location


                flash("Found an A1 region in %s" % record.name)
                location = re.search(r"\d*:\d*", str(best_location))

                setattr(record, "a1", best_seq)
                setattr(record, "a1_loc", location[0])
                db.session.add(record)
                db.session.commit()

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')

    @action('b_check_a2', 'Check for A2 region')
    def action_check_a2(self, ids):
        try:
            best_align = 0
            query = SequenceRecords.query.filter(SequenceRecords.uid.in_(ids))
            for record in query.all():
                print(record.name)
                seq_record = bio_db.lookup(primary_id=record.name)

                for feature in seq_record.features:
                    if 'translation' in feature.qualifiers:
                        alignment = pairwise2.align.globalms(yenA2, feature.qualifiers['translation'][0],
                                                             2, -1, -.5, -.1, score_only=True)
                        if alignment > best_align:
                            best_align = alignment
                            best_seq = feature.qualifiers['translation'][0]
                            best_location = feature.location


                flash("Found an A2 region in %s" % record.name)
                location = re.search(r"\d*:\d*", str(best_location))
                setattr(record, "a2", best_seq)
                setattr(record, "a2_loc", location[0])
                db.session.add(record)
                db.session.commit()

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')

    @action('d_buildprofile_a2', 'Build profile based on A2')
    def action_build_profile_a2(self, ids):
        try:
            query = SequenceRecords.query.filter(SequenceRecords.uid.in_(ids))
            align_list = []
            for record in query.all():
                if record.a2 == "Not tested":
                    pass
                else:
                    align_record = SeqRecord(Seq(record.a2, generic_protein), id=str(record.name) + "_" + "A2")
                    align_list.append(align_record)

                    # Write sequences to FASTA

            SeqIO.write(align_list, "align.fasta", "fasta")
            muscle_cline = MuscleCommandline(input="align.fasta")
            child = subprocess.Popen(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     universal_newlines=True, shell=(sys.platform != "win32"))
            alignment = AlignIO.read(child.stdout, "fasta")
            AlignIO.write(alignment, "align.aln", "fasta")
            file = subprocess.call(["hmmbuild", "profile.hmm", "align.aln"])

            file_contents = file.stream.read().decode("utf-8")

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')


class ProfileView(ModelView):

    create_modal = True
    edit_modal = True
    can_create = False
    can_view_details = True



# Form for uploading files
class UploadForm(FlaskForm):
    name = StringField('What ID should we give this feature?', validators=[DataRequired("Not completed")])
    file = FileField('Upload the FASTA file', validators=[DataRequired("Not completed")])

    upload_submit = SubmitField("Upload files")

# View for uploading files
class UploadView(BaseView):

    @expose("/", methods =('GET', 'POST'))
    def upload(self):
        form = UploadForm()

        if request.method == 'POST':
            # Get the file and ID to name these sequences
            region = form.name.data
            filename = allfiles.save(request.files['file'])
            # Create the initial seqRecords
            phyloisland.seqDict = {}
            phyloisland.unmappable = []
            phyloisland.defaultValue("static/uploads/" + filename, region)

            records = phyloisland.seqDict

            for record in records:
                current = records[record]
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

                entry = SequenceRecords(name, species, strain, description, a1, a1_loc, a2, a2_loc, overlap, distance, sequence)
                check = SequenceRecords.query.filter_by(name=name).first()

                if check:
                    if region in current.annotations.keys():
                        print("YEAH BABY")
                    else:
                        print("NAH BOY")
                    setattr(check, region.lower(), current.annotations[region] if region in current.annotations.keys() else "Not tested")
                    setattr(check, region.lower() + "_loc", current.annotations[region + "_location"] if region + "_location" in current.annotations.keys() else "Not tested")
                    db.session.add(check)

                else:
                    seq_list = []
                    db.session.add(entry)
                    seq_list.append(current)
                    bio_db.load(seq_list)
                    bio_server.commit()
                db.session.commit()

            return self.render("upload_admin.html", form=form, records=records)
        return self.render("upload_admin.html", form=form)

class MyHomeView(AdminIndexView):
	@expose('/phyloisland_experimental')
	def index(self):
		return self.render('admin/index.html')


@application.route(local('/'))
def index():
    return '<a href="/admin/">Click me to get to Phylo Island please!</a>'



#admin = Admin(application,name="Phylo Island", template_mode="bootstrap3")
admin = Admin(application, index_view=AdminIndexView(name='Experimental', url="/phyloisland_experimental"))
#admin = Admin(application, index_view=MyHomeView())
admin.add_view(UploadView(name='Upload', endpoint='upload_admin'))
admin.add_view(SequenceRecordsView(SequenceRecords, db.session, endpoint="seq_view"))  # working version
admin.add_view(ProfileView(Profiles, db.session, endpoint="profiles"))
if __name__ == "__main__":
	application.run(debug=True, host='0.0.0.0')

