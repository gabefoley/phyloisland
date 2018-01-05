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
import servers


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



a1_profile_path = "tmp/A1_profile.hmm"
a2_profile_path = "tmp/A2_profile.hmm"



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


    def apply(self, query, value, alias="None"):

        species_list = []
        id_list = []

        for record in servers.bio_db.values():
            species = (" ".join(record.annotations.get('organism').split()[0:2]))
            if species in species_list:
                continue
            else:
                species_list.append(species)
                id_list.append(record.id)
        return query.filter(self.get_column(alias).in_(id_list))

    def operation(self):
        return 'Yes'


class SequenceRecordsView(ModelView):
    create_modal = True
    edit_modal = True
    can_create = False
    can_view_details = True

    def _seqdescription_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window


        if model.sequence:
            return model.sequence[:15] + "..."
        else:
            return model.sequence
    column_formatters = {
        'sequence': _seqdescription_formatter,
    }

    @action('set_A1_reference', 'Set this sequence as the A1 reference')
    def action_set_A1_reference(self, ids):
        if len(ids) > 1:
            flash('Only select a single record')
        else:
            query = models.SequenceRecords.query.filter(models.SequenceRecords.uid.in_(ids))
            for record in query.all():

                # Check for a previous reference sequence
                old_genome_reference = models.GenomeRecords.query.filter_by(a1_ref=1).first()
                old_seq_reference = models.SequenceRecords.query.filter_by(a1_ref=1).first()

                if (old_genome_reference):
                    # Remove the previous reference sequence
                    setattr(old_genome_reference, "a1_ref", 0)
                    servers.db.session.add(old_genome_reference)

                if (old_seq_reference):
                    # Remove the previous reference sequence
                    setattr(old_seq_reference, "a1_ref", 0)
                    servers.db.session.add(old_seq_reference)

                # Set the new reference sequence
                setattr(record, "a1_ref", 1)

                # Commit the changed record
                servers.db.session.add(record)
                servers.db.session.commit()

                # global a1_reference
                # a1_reference = SeqRecord(Seq(record.sequence, generic_protein), id=str(record.name) + "_reference")
                flash ("The A1 region from %s has been set as the reference A1 sequence" % record.name)


    @action('set_A2_reference', 'Set this sequence as the A2 reference')
    def action_set_A2_reference(self, ids):
        if len(ids) > 1:
            flash('Only select a single record')
        else:
            query = models.SequenceRecords.query.filter(models.SequenceRecords.uid.in_(ids))
            for record in query.all():

                # Check for a previous reference sequence
                old_genome_reference = models.GenomeRecords.query.filter_by(a2_ref=1).first()
                old_seq_reference = models.SequenceRecords.query.filter_by(a2_ref=1).first()

                if (old_genome_reference):
                    # Remove the previous reference sequence
                    setattr(old_genome_reference, "a2_ref", 0)
                    servers.db.session.add(old_genome_reference)

                if (old_seq_reference):
                    # Remove the previous reference sequence
                    setattr(old_seq_reference, "a2_ref", 0)
                    servers.db.session.add(old_seq_reference)

                # Set the new reference sequence
                setattr(record, "a2_ref", 1)

                # Commit the changed record
                servers.db.session.add(record)
                servers.db.session.commit()

                flash ("The A2 region from %s has been set as the reference A2 sequence" % record.name)



    @action('generate_profile', 'Generate a profile from these sequences')
    def action_generate_profile(self, ids):
        try:
            query = models.SequenceRecords.query.filter(models.SequenceRecords.uid.in_(ids))
            align_list = []
            for record in query.all():
                align_record = SeqRecord(Seq(record.sequence, generic_protein), id=str(record.name) + "_" + "seq")
                align_list.append(align_record)

            createProfile(align_list)


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext(ex))


class GenomeRecordsView(ModelView):
    column_list = ('name', 'species', 'strain', 'description', 'a1_ref', 'a2_ref', 'a1', 'a1_length', 'a1_loc', 'a2', 'a2_length', 'a2_loc', 'overlap', 'distance', 'sequence')
    column_searchable_list = ['name', 'species', 'a1', 'a2', 'overlap']
    create_modal = True
    edit_modal = True
    can_create = False
    can_view_details = True

    def _a1description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.a1:
            return model.a1[:15] + "..."
        else:
            return model.a1
    def _a2description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.a2:
            return model.a2[:15] + "..."
        else:
            return model.a2


    def _seqdescription_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window


        if model.sequence:
            return model.sequence[:15] + "..."
        else:
            return model.sequence
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
                          models.GenomeRecords.name, 'Get unique species'),
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



    @action('item5_getoverlap', 'Get Overlap')
    def item5_getoverlap(self, ids):
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

            flash(gettext('Failed to get overlap %(error)s', error=str(ex)), 'error')

    @action('item6_A1_reference', 'Set this A1 region as the reference region')
    def item6_set_a1_reference(self, ids):
        if len(ids) > 1:
            flash('Only select a single record')
        else:
            query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
            for record in query.all():
                if (record.a1):

                    # Check for a previous reference sequence
                    old_genome_reference = models.GenomeRecords.query.filter_by(a1_ref=1).first()
                    old_seq_reference = models.SequenceRecords.query.filter_by(a1_ref=1).first()


                    if (old_genome_reference):
                        # Remove the previous reference sequence
                        setattr(old_genome_reference, "a1_ref", 0)
                        servers.db.session.add(old_genome_reference)

                    if (old_seq_reference):
                        # Remove the previous reference sequence
                        setattr(old_seq_reference, "a1_ref", 0)
                        servers.db.session.add(old_seq_reference)

                    # Set the new reference sequence
                    setattr(record, "a1_ref", 1)

                    # Commit the changed record
                    servers.db.session.add(record)
                    servers.db.session.commit()
                    flash ("The A1 region from %s has been set as the reference A1 sequence" % record.name)
                else:
                    flash ("That record doesn't have an A1 region associated with it")


    @action('item7_set_A2_reference', 'Set this A2 region as the reference region')
    def item7_set_a2_reference(self, ids):
        if len(ids) > 1:
            flash('Only select a single record')
        else:
            query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
            for record in query.all():
                if (record.a2):

                    # Check for a previous reference sequence
                    old_genome_reference = models.GenomeRecords.query.filter_by(a2_ref=1).first()
                    old_seq_reference = models.SequenceRecords.query.filter_by(a2_ref=1).first()

                    if (old_genome_reference):
                        # Remove the previous reference sequence
                        setattr(old_genome_reference, "a2_ref", 0)
                        servers.db.session.add(old_genome_reference)

                    if (old_seq_reference):
                        # Remove the previous reference sequence
                        setattr(old_seq_reference, "a2_ref", 0)
                        servers.db.session.add(old_seq_reference)

                    # Set the new reference sequence
                    setattr(record, "a2_ref", 1)

                    # Commit the changed record
                    servers.db.session.add(record)
                    servers.db.session.commit()
                    flash ("The A2 region from %s has been set as the reference A2 sequence" % record.name)
                else:
                    flash ("That record doesn't have an A2 region associated with it")

    @action('item_a1_profile_align_a1', 'Check for A1 region with a profile')
    def item_a1_profile_align_a1(self, ids):
        try:

            # Check if a reference profile for A1 exists
            a1_profile_reference = models.Profile.query.filter_by(a1_profile_ref=1).first()

            print ('looking for alibrandi / a1')
            print (a1_profile_reference)

            if (a1_profile_reference):
                print ("Using the A1 profile named %s to check for A1 regions" % (a1_profile_reference.name))
                checkForFeature.get_feature_location_with_profile(ids, a1_profile_path, "a1", "a1_loc")
            else:
                flash("Please set a profile as the A1 reference profile first", "error")

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for A1 region -  %(error)s', error=str(ex)), 'error')


    @action('item_a2_profile_align_a2', 'Check for A2 region with a profile')
    def item_a2_profile_align_a2(self, ids):
        try:

            # Check if a reference profile for A1 exists
            a2_profile_reference = models.Profile.query.filter_by(a2_profile_ref=1).first()

            if (a2_profile_reference):
                print ("Using the A2 profile named %s to check for A2 regions" % (a2_profile_reference.name))
                checkForFeature.get_feature_location_with_profile(ids, a2_profile_path, "a2", "a2_loc")
            else:
                flash("Please set a profile as the A2 reference profile first", "error")

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for A1 region -  %(error)s', error=str(ex)), 'error')

    @action('item3_check_a1', 'Check for A1 region')
    def item3_check_a1(self, ids):
        try:

            checkForFeature.getFeatureLocation(ids, a1_reference, "a1", "a1_loc", "a1_length")

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for A1 region -  %(error)s', error=str(ex)), 'error')

    @action('item4_check_a2', 'Check for A2 region')
    def item4_check_a2(self, ids):
        try:

            checkForFeature.getFeatureLocation(ids, a2_reference, "a2", "a2_loc", "a2_length")




        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for A2 region -  %(error)s', error=str(ex)), 'error')



    @action('item1_delete_a1', 'Delete A1 region')
    def item1_delete_a1(self, ids):
        try:

            checkForFeature.deleteFeature(ids, "a1", "a1_loc", "a1_length")




        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to delete region. %(error)s', error=str(ex)), 'error')


    @action('item2_delete_a2', 'Delete A2 region')
    def item2_delete_a2(self, ids):
        try:

            checkForFeature.deleteFeature(ids, "a2", "a2_loc", "a2_length")



        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to delete region. %(error)s', error=str(ex)), 'error')

    @action('item8_generate_profile_a2', 'Generate profile based on A1')
    def item8_generate_profile_a1(self, ids):
        try:
            query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
            align_list = []
            for record in query.all():
                if record.a1 == "":
                    pass
                else:
                    align_record = SeqRecord(Seq(record.a1, generic_protein), id=str(record.name) + "_" + "A1")
                    align_list.append(align_record)

            createProfile(align_list)


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to generate profile based on A1. %(error)s', error=str(ex)), 'error')

    @action('item9_generate_profile_a2', 'Generate profile based on A2')
    def item9_generate_profile_a2(self, ids):
        try:
            query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
            align_list = []
            for record in query.all():
                if record.a2 == "":
                    pass
                else:
                    align_record = SeqRecord(Seq(record.a2, generic_protein), id=str(record.name) + "_" + "A2")
                    align_list.append(align_record)

            createProfile(align_list)


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to generate profile based on A2. %(error)s', error=str(ex)), 'error')


class ProfileView(ModelView):
    """
    View of the Profile database for storing HMM Profiles """
    column_list = ('name', 'a1_profile_ref', 'a2_profile_ref', 'download')
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

    @action('item1_set_A1_reference', 'Set this profile as the A1 reference profile')
    def item1_set_a1_reference(self, ids):
        if len(ids) > 1:
            flash('Only select a single record')
        else:
            query = models.Profile.query.filter(models.Profile.uid.in_(ids))
            for record in query.all():

                # Check for a previous reference profile
                old_profile_reference = models.Profile.query.filter_by(a1_profile_ref=1).first()

                if (old_profile_reference):
                    # Remove the previous reference profile
                    setattr(old_profile_reference, "a1_profile_ref", 0)
                    servers.db.session.add(old_profile_reference)


                # Set the new reference profile
                setattr(record, "a1_profile_ref", 1)

                # Commit the changed record
                servers.db.session.add(record)
                servers.db.session.commit()

                # Write the new profile to the tmp folder ready to be used




                with open(a1_profile_path, 'w') as profile_path:
                    print (record)
                    # print (record.profile)
                    profile_path.write(record.profile.decode('utf-8'))

                # global a1_profile_reference
                #
                # a1_profile_reference["tmp/A1_profile.hmm"] = record.name
                flash ("The profile named %s has been set as the reference A1 profile" % record.name)

    @action('item2_set_A1_reference', 'Set this profile as the A2 reference profile')
    def item2_set_a1_reference(self, ids):
        if len(ids) > 1:
            flash('Only select a single record')
        else:
            query = models.Profile.query.filter(models.Profile.uid.in_(ids))
            for record in query.all():
                for record in query.all():

                    # Check for a previous reference profile
                    old_profile_reference = models.Profile.query.filter_by(a2_profile_ref=1).first()

                    if (old_profile_reference):
                        # Remove the previous reference profile
                        setattr(old_profile_reference, "a2_profile_ref", 0)
                        servers.db.session.add(old_profile_reference)

                    # Set the new reference profile
                    setattr(record, "a2_profile_ref", 1)

                    # Commit the changed record
                    servers.db.session.add(record)
                    servers.db.session.commit()

                    # Write the new profile to the tmp folder ready to be used
                    with open(a2_profile_path, 'w') as profile_path:
                        print(record)
                        print(record.profile)
                        profile_path.write(record.profile.decode('utf-8'))

                    flash("The profile named %s has been set as the reference A2 profile" % record.name)



class UploadView(BaseView):
    """
    View for uploading files
    """

    @expose("/", methods =('GET', 'POST'))
    def upload(self):
        form = models.UploadForm()

        if request.method == 'POST' and form.validate():


            # Get the sequences
            filename = servers.allfiles.save(request.files['file'])

            # Create the initial seqRecords
            phyloisland.seqDict = {}
            phyloisland.unmappable = []

            seq_records = SeqIO.parse("static/uploads/" + filename, "fasta")

            for record in seq_records:
                seq_name = record.id
                seq_description = record.description.split(">")[0]
                seq_species = seq_description.split("[")[1].split("]")[0]
                seq_sequence = str(record.seq)

                seq_entry = models.SequenceRecords(seq_name, seq_species, seq_description, seq_sequence, 0, 0)

                # Check if the sequence record already exists
                seq_check = models.SequenceRecords.query.filter_by(name=seq_name).first()
                if seq_check == None:
                    print ('Adding sequence with ID - %s from species - %s to the sequence database' % (seq_name, seq_species) + "\n")
                    servers.db.session.add(seq_entry)
                    servers.db.session.commit()
                else:
                    print ('Sequence with ID - %s from species - %s already exists in the sequence database' % (seq_name, seq_species) + "\n")




            # Map the sequence to its genome
            genomeResults = mapToGenome.getFullGenome("static/uploads/" + filename)

            for genome in genomeResults:
                print ('Could not find')
                print (genome)

                flash("Couldn't find " + str(genome))

            records = mapToGenome.seqDict

            for record in records:


                current = records[record]
                name = current.id



                species = " ".join(current.annotations.get('organism').split()[0:2])
                strain = current.annotations['source']
                sequence = str(current.seq)
                print (current)
                print (species)

                description = current.description
                a1 = current.annotations["A1"] if "A1" in current.annotations.keys() else ""
                a1_length = None
                a1_loc = current.annotations["A1_location"] if "A1_location" in current.annotations.keys() else ""
                a2 = current.annotations["A2"] if "A2" in current.annotations.keys() else ""
                a2_length = None
                a2_loc = current.annotations["A2_location"] if "A2_location" in current.annotations.keys() else ""
                overlap = current.annotations["Overlap"] if "Overlap" in current.annotations.keys() else ""
                distance = ""

                entry = models.GenomeRecords(name, species, strain, description, a1, a1_length, a1_loc, a2, a2_length, a2_loc, overlap, distance, sequence, 0, 0)
                check = models.GenomeRecords.query.filter_by(name=name).first()

                # Check to see if the genome record already exists
                if check:
                    print ("The genome record - %s from species - %s already exists in the database" % (name, species))

                    continue
                #     # if region in current.annotations.keys():
                #     #     continue
                #     # else:
                #     #     setattr(check, region.lower(), current.annotations[region] if region in current.annotations.keys() else "")
                #     #     setattr(check, region.lower() + "_loc", current.annotations[region + "_location"] if region + "_location" in current.annotations.keys() else "")
                #     #     servers.db.session.add(check)

                else:
                    print ("Adding the genome record - %s from species - %s to the genome database" % (name, species))

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

def createProfile(align_list):
    SeqIO.write(align_list, "align.fasta", "fasta")
    muscle_cline = MuscleCommandline(input="align.fasta")
    # result = subprocess.run(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True) )
    child = subprocess.Popen(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True, shell=(sys.platform != "win32"))
    child.wait()

    alignment = AlignIO.read(child.stdout, "fasta")
    AlignIO.write(alignment, "align.aln", "fasta")
    result = subprocess.call(["hmmbuild", "profile3.hmm", "align.aln"], stdout=subprocess.PIPE)
    file = open('profile3.hmm', 'rb')

    saveProfile(file)



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
