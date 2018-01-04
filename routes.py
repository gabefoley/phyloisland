from os.path import join
from flask import flash
from flask_admin import Admin, BaseView, expose
from flask_admin.contrib.sqla import ModelView
import os
import phyloisland
from flask_admin.actions import action
import gettext
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
import subprocess
import sys
from Bio.SeqRecord import SeqRecord
from flask import request
from flask_admin.contrib.sqla import filters
from flask_admin.contrib.sqla.filters import BaseSQLAFilter
from flask_admin import AdminIndexView
import os.path as op
import mapToGenome
import checkForFeature
import models, servers
import io
from flask import send_file
from markupsafe import Markup



try:
    from wtforms.fields.core import _unset_value as unset_value
except ImportError:
    from wtforms.utils import unset_value


# Create directory for file fields to use
file_path = op.join(op.dirname(__file__), 'filesdir')
try:
    os.mkdir(file_path)
except OSError:
    pass
BASE_ROUTE = '/' + servers.base_route

def local(route: str) -> str:
    if BASE_ROUTE == '/':
        return route
    else:
        return join(BASE_ROUTE, route[1:])


a1_reference = SeqRecord(Seq("MDKYNNYSNVIKNKSSISPLLAAAAKIEPEITVLSSASKSNRSQYSQSLADTLLGLGYRSIFDIAKVSRQRFIKRHDESLLGNGAVIFDKAVSMANQVLQKYRKNRLEKSNSPLVPQTSSSTDASSESQTNKLPEYNQLFPEPWDNFCRPGAIEALDSPASYLLDLYKFIQSVELDGSNQARKLETRRADIPKLSLDNDALYKEVTALSIVNDVLSGSAREYIDQSGQADKAVNQILGDTHFPFTLPYSLPTQQINKGLGASNIELGTVIQRVDPQFSWNTTQEKYNQVLLAYTQLSSEQIALLSLPDVFTQNFLTQTELSAGYLSASTTEILAEKDLSRHGYIVKAADNIKGPTQLVEHSDASYDVIELTCTNQAKETITVKLRGENIITYQRTKARMVPFDNSSPFSRQLKLTFVAEDNPSLGNLDKGPYFANMDIYAAEWVRENVSSETMVSRPFLTMTYRIAIAKAGASLEELQPEADAFFINNFGLSAEDSSQLVKLVAFGDQTGSKAEEIESLLSCGENLPIVSPNVIFANPIFGSYFNDEPFPAPYHFGGVYINAHQRNAMTIIRAEGGREIQSLSNFRLERLNRFIRLQRWLDLPSHQLDLLLTSVMQADADNSQQEITEPVLKSLGLFRHLNLQYKITPEIFSSWLYQLTPFAVSGEIAFFDRIFNREQLFDQPFILDGGSFTYLDAKGSDAKSVKQLCAGLNISAVTFQFIAPLVQSALGLEAGTLVRSFEVVSSLYRLVSIPQTFGLSTEDGLILMNILTDEMGYLAKQPAFDDKQTQDKDFLSIILKMEALSAWLTKNNLTPASLALLLGVTRLAVVPTNNMVTFFKGIANGLSENVCLTTDDFQRQELEGADWWTLLSTNQVIDDMGLVLDIHPVWGKSDEEMLMEKIQSIGVSNDNNTLSIIVQILIQAKNAQENLLSQTISAEYGVERSVVPLQLRWLGSNVYSVLNQVLNNTPTDISSIVPKLSELTYSLLIYTQLINSLKLNKEFIFLRLTQPNWLGLTQPKLSTQLSLPEIYLITCYQDWVVNANKNEDSIHEYLEFANIKKTEAEKTLVDNSEKCAELLAEILAWDAGEILKAASLLGLNPPQATNVFEIDWIRRLQTLSEKTMISTEYLWQMGDLTENSEFSLKEGVGEAVMAALKAQGDSDNV", generic_protein), id="Yersinia entomophaga A1")

a2_reference = SeqRecord(Seq("MSNSIEAKLQEDLRDALVDYYLGQIVPNSKDFTNLRSTIKNVDDLYDHLLLDTQVSAKVITSRLSLVTQSVQQYINRIALNLEPGLSINQQEATDWEEFANRYGYWAANQQLRMFPEIYVDPTLRLTKTEFFFQLESALNQGKLTDDVAQKAVLGYLNNFEEVSNLEIIAGYQDGIDIENDKTYFVARTRMQPYRYFWRSLDASQRNANSQELYPTAWSEWKAISVPLENVANGIVRPIMMDNRLYISWFEVAEEKETDSDGNIIVSGRYRTKIRLAHLGFDGVWSSGTTLREEVLADQMEEMIAVVDRMEDEPRLALVAFKEMSESWDVVFSYICDSMLIESSNLPTTTHPPKPGDGDKGLSDLDDYGANLVWFYLHETANGGKAEYKQLILYPVIINRDWPIELDKTHQGDFGTVDDFTLNSNYTGDELSLYLQSSSTYKYDFSKSKNIIYGIWKEDANNNRCWLNYKLLTPEDYEPQINATLVMCDKGDVNIITGFSLPNGGVDAGGKIKVTLRVGKKLRDKFQIKQFSQTQYLQFPEASSADVWYIGKQIRLNTLFAKELIGKASRSLDLVLSWETQNSRLEEAILGGAAELIDLDGANGIYFWELFFHMPFMVSWRFNVEQRYEDANRWVKYLFNPFECEDEPALLLGKPPYWNSRPLVDEPFKGYSLTQPSDPDAIAASDPIHYRKAVFNFLTKNIIDQGDMEYRKLQPSARTLARLSYSTASSLLGRRPDVQLTSFWQPLTLEDASYKTDSEIRAIEMQSQPLTFEPVVHDQTMSAVDNDIFMYPMNNELRGLWDRIENRIYNLRHNLTLDGKEINMDLYDSSISPRGLMKQRYQRVVTARNASKMNFKVPNYRFEPMLNRSKSGVETLIQFGSTLLSLLERKDSLSFDAYQMIQSGDLYRFSIDLQQQDIDINKASLEALQVSKQSAQDRYDHFKELYDENISSTEQKVIELQSQAANSLLMAQGMRTAAAALDVIPNIYGLAVGGSHWGAPLNAAAEIIMIKYQADSSKSESLSVSESYRRRRQEWELQYKQAEWEVNSVEQQINLQNMQIKAANKRLEQVEAQQQQAMALLDYFSERFTNESLYTWLISQLSSLYLQAYDAVLSLCLSAEASLLYELNLGEQSFVGGGGWNDLYQGLMAGETLKLALMRMERVYVEQNSRRQEITKTISLKALLGESWPAELNKLKQKTPINFNLEEQIFVEDYQELYQRRIKSVSVSLPMLVGPYEDVCAQLTQTSSSYSTRADLKTVENMLTKRTFADTPHLVRSIQPNQQISLSTGVNDSGLFMLNFDDERFLPFEGSGVDSSWRLQFTNLKQNLDSLNDVILHVKYTAAIGSSTFSQGVRKILANINNDE", generic_protein), id="Yersinia entomophaga A2")


class ProfileView(ModelView):
    """
    View of the Profile database for storing HMM Profiles """
    column_list = ('name', 'size', 'filename', 'mimetype', 'download')
    form_columns = ('name', 'profile')

    form_extra_fields = {'profile': models.BlobUploadField(
        label='File',
        allowed_extensions=['pdf', 'doc', 'docx', 'xls', 'xlsx', 'png', 'jpg', 'jpeg', 'gif', 'hmm', 'fa', 'fasta'],
        size_field='size',
        filename_field='filename',
        mimetype_field='mimetype'
    )}

    def _download_formatter(self, context, model, name):
        return Markup(
            "<a href='{url}' target='_blank'>Download</a>".format(url=self.get_url('download_blob', id=model.uid)))

    column_formatters = {
        'download': _download_formatter,
    }


@servers.application.route("/" + servers.base_route + "/<int:id>", methods=['GET'])
def download_blob(id):
    """
    Route for downloading profiles from the Profile view
    :param id: Profile to download
    :return:
    """
    _profile = models.Profile.query.get_or_404(id)
    return send_file(
        io.BytesIO(_profile.profile),
        attachment_filename=_profile.filename,
        mimetype=_profile.mimetype
    )


class FilterInAListMaybe(BaseSQLAFilter):

    def apply(self, query, value, alias=None):
        ids = {}
        unique_ids = []
        unique_items = []

        for id in ids:
            if id in unique_ids:
                continue
            else:
                unique_ids.add(id)
                unique_items.add()


        listup = ['Bradyrhizobium diazoefficiens', 'FMWJ01000001.1', 'CPYD01000004.1', 'DQ400808.1']
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
    create_modal = True
    edit_modal = True
    can_create = False
    can_view_details = True

    @action('set_A1_reference', 'Set this sequence as the A1 reference')
    def action_set_A1_reference(self, ids):
        if len(ids) > 1:
            flash('Only select a single record')
        else:
            query = models.SequenceRecords.query.filter(models.SequenceRecords.uid.in_(ids))
            for record in query.all():
                global a1_reference
                a1_reference = SeqRecord(Seq(record.sequence, generic_protein), id=str(record.name) + "_reference")
                flash ("The A1 region from %s has been set as the reference A1 sequence" % record.name)


    @action('set_A2_reference', 'Set this sequence as the A2 reference')
    def action_set_A2_reference(self, ids):
        if len(ids) > 1:
            flash('Only select a single record')
        else:
            query = models.SequenceRecords.query.filter(models.SequenceRecords.uid.in_(ids))
            for record in query.all():
                global a2_reference
                a2_reference = SeqRecord(Seq(record.sequence, generic_protein), id=str(record.name) + "_reference")
                flash ("The A1 region from %s has been set as the reference A2 sequence" % record.name)

class GenomeRecordsView(ModelView):
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
                      filters.FilterLike(models.GenomeRecords.overlap, 'Overlap',
                                         options=(('True', 'True'), ('False', 'False'))),
                      'sequence',
                      FilterInAListMaybe(
                          models.GenomeRecords.species, 'Get unique species'),
                      filters.IntGreaterFilter(models.GenomeRecords.distance, 'Distance greater than'),
                      filters.IntSmallerFilter(models.GenomeRecords.distance, 'Distance smaller than')

                      )

    def after_model_delete(self, model):
        """
        Function to delete the record from the BioSQL database as well
        :param model: Record to delete
        :return:
        """
        del servers.bio_db[model.uid]
        servers.bio_server.commit()



    @action('c_getoverlap', 'Get Overlap')
    def action_getoverlap(self, ids):
        try:
            query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
            for record in query.all():
                if record.a1 == "" or record.a2 == "":
                    pass
                else:
                    distance = str(phyloisland.getDistance(record.a1_loc, record.a2_loc))
                    print (distance)
                    if len(distance) > 1:
                        record.overlap = "False"
                        record.distance = distance
                    else:
                        record.overlap = "True"
                        record.distance = ""

                    servers.db.session.add(record)
                    servers.db.session.commit()

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')

    @action('set_A1_reference', 'Set this A1 region as the reference region')
    def action_set_A1_reference(self, ids):
        if len(ids) > 1:
            flash('Only select a single record')
        else:
            query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
            for record in query.all():
                if (record.a2):
                    global a1_reference
                    a1_reference = SeqRecord(Seq(record.a1, generic_protein), id=str(record.name) + "_reference")
                    flash ("The A1 region from %s has been set as the reference A1 sequence" % record.name)
                else:
                    flash ("That record doesn't have an A1 region associated with it")


    @action('set_A2_reference', 'Set this A2 region as the reference region')
    def action_set_A2_reference(self, ids):
        if len(ids) > 1:
            flash('Only select a single record')
        else:
            query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
            for record in query.all():
                if (record.a2):
                    global a2_reference
                    a2_reference = SeqRecord(Seq(record.a2, generic_protein), id=str(record.name) + "_reference")
                    flash ("The A2 region from %s has been set as the reference A2 sequence" % record.name)
                else:
                    flash ("That record doesn't have an A2 region associated with it")


    @action('profile_align', 'Check for A2 region with a profile')
    def action_profile_align_a2(self, ids):
        try:

            checkForFeature.get_feature_location_with_profile(ids, "refs/testprofile.hmm", "a2", "a2_loc")

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')


    @action('a_check_a1', 'Check for A1 region')
    def action_check_a1(self, ids):
        try:

            checkForFeature.getFeatureLocation(ids, a1_reference, "a1", "a1_loc")

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')

    @action('b_check_a2', 'Check for A2 region')
    def action_check_a2(self, ids):
        try:

            #TODO: This obviously shouldn't be hard coded either!
            checkForFeature.getFeatureLocation(ids, a2_reference, "a2", "a2_loc")




        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to approve users. %(error)s', error=str(ex)), 'error')



    @action('c_delete_a1', 'Delete A1 region')
    def action_delete_a1(self, ids):
        try:

            checkForFeature.deleteFeature(ids, "a1", "a1_loc")




        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to delete records. %(error)s', error=str(ex)), 'error')


    @action('c_delete_a2', 'Delete A2 region')
    def action_delete_a2(self, ids):
        try:

            checkForFeature.deleteFeature(ids, "a2", "a2_loc")



        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to delete records. %(error)s', error=str(ex)), 'error')




    @action('d_buildprofile_a2', 'Build profile based on A2')
    def action_build_profile_a2(self, ids):
        try:
            query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
            align_list = []
            for record in query.all():
                if record.a2 == "":
                    pass
                else:
                    align_record = SeqRecord(Seq(record.a2, generic_protein), id=str(record.name) + "_" + "A2")
                    align_list.append(align_record)

                    # Write sequences to FASTA

            SeqIO.write(align_list, "align.fasta", "fasta")
            muscle_cline = MuscleCommandline(input="align.fasta")

            print (type(muscle_cline))
            # result = subprocess.run(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True) )
            child = subprocess.Popen(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     universal_newlines=True, shell=(sys.platform != "win32"))
            child.wait()


            alignment = AlignIO.read(child.stdout, "fasta")
            AlignIO.write(alignment, "align.aln", "fasta")
            result = subprocess.call(["hmmbuild", "profile3.hmm", "align.aln"], stdout=subprocess.PIPE)
            file = open('profile3.hmm', 'rb')


            saveProfile(file)

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext(ex))




class UploadView(BaseView):
    """
    View for uploading files
    """

    @expose("/", methods =('GET', 'POST'))
    def upload(self):
        form = models.UploadForm()

        if request.method == 'POST':

            # Get the sequences
            filename = servers.allfiles.save(request.files['file'])

            # Create the initial seqRecords
            phyloisland.seqDict = {}
            phyloisland.unmappable = []

            seq_records = SeqIO.parse("static/uploads/" + filename, "fasta")

            print (seq_records)

            for record in seq_records:
                print (record)
                print (record.annotations)


                seq_name = record.id
                seq_description = record.description.split(">")[0]
                seq_species = seq_description.split("[")[1].split("]")[0]
                seq_sequence = str(record.seq)

                seq_entry = models.SequenceRecords(seq_name, seq_species, seq_description, seq_sequence)

                # Check if the sequence record already exists
                seq_check = models.GenomeRecords.query.filter_by(name=seq_name).first()
                if not seq_check:
                    servers.db.session.add(seq_entry)
                    servers.db.session.commit()


                print (seq_name)
                print (seq_description)
                print (seq_species)

            genomeResults = mapToGenome.getFullGenome("static/uploads/" + filename)

            print (genomeResults)

            for genome in genomeResults:
                # print (genome)
                # print (type(genome))
                flash("Couldn't find " + str(genome))
                # flash(gettext("Couldn't find").join(str(genome)), 'error')


            # flash(gettext("Couldn't find").join([v for v in genomeResults]), 'error')

            # phyloisland.defaultValue("static/uploads/" + filename, region)

            records = mapToGenome.seqDict

            for record in records:


                current = records[record]
                name = current.id

                species = current.annotations.get('organism')
                strain = current.annotations['source']
                sequence = str(current.seq)

                description = current.description
                a1 = current.annotations["A1"] if "A1" in current.annotations.keys() else ""
                a1_loc = current.annotations["A1_location"] if "A1_location" in current.annotations.keys() else ""
                a2 = current.annotations["A2"] if "A2" in current.annotations.keys() else ""
                a2_loc = current.annotations["A2_location"] if "A2_location" in current.annotations.keys() else ""
                overlap = current.annotations["Overlap"] if "Overlap" in current.annotations.keys() else ""
                distance = ""

                print ('and the record id is ')
                print (name)
                print ('and the seq name is ')
                print (seq_name)


                entry = models.GenomeRecords(name, species, strain, description, a1, a1_loc, a2, a2_loc, overlap, distance, sequence)
                check = models.GenomeRecords.query.filter_by(name=name).first()
                print ('lets check')
                print (check)
                # Check to see if the genome record already exists
                if check:
                    continue
                #     # if region in current.annotations.keys():
                #     #     continue
                #     # else:
                #     #     setattr(check, region.lower(), current.annotations[region] if region in current.annotations.keys() else "")
                #     #     setattr(check, region.lower() + "_loc", current.annotations[region + "_location"] if region + "_location" in current.annotations.keys() else "")
                #     #     servers.db.session.add(check)

                else:
                    seq_list = []
                    servers.db.session.add(entry)
                    seq_list.append(current)
                    servers.bio_db.load(seq_list)
                    servers.bio_server.commit()
                    servers.db.session.commit()



            return self.render("upload_admin.html", form=form, records=records)
        return self.render("upload_admin.html", form=form)

class MyHomeView(AdminIndexView):
	@expose(servers.base_route)
	def index(self):
		return self.render('admin/index.html')


def saveProfile(profile):
    """
    Save a profile into the database
    :param profile: Profile to save
    :return:
    """

    name = phyloisland.randstring(5)
    blobMix = models.BlobMixin("application/octet-stream", name, profile.read(), '566666')
    profileEntry = models.Profile(name)
    profileEntry.set_blobMix(blobMix)
    servers.db.session.add(profileEntry)
    servers.db.session.commit()


# Setup the main flask-admin application
admin = Admin(servers.application, index_view=AdminIndexView(name= servers.base_name, url="/" + servers.base_route), template_mode="bootstrap3")

# Add the views to the flask-admin application
admin.add_view(UploadView(name='Upload', endpoint='upload_admin'))
admin.add_view(SequenceRecordsView(models.SequenceRecords, servers.db.session, endpoint="seq_view"))  # working version

admin.add_view(GenomeRecordsView(models.GenomeRecords, servers.db.session, endpoint="genome_view"))  # working version
admin.add_view(ProfileView(model=models.Profile, session=servers.db.session, name='Profiles'))

if __name__ == '__main__':
    servers.application.run(debug=True, port=7777)
