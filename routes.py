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
import utilities
import time

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

genome_records = []


def setReferenceProfile(region):
    # filter = globals()[region + '_profile_ref=1']()
    query = eval('models.Profile.query.filter_by(' + region + '_profile_ref=1).first()')
    if query:
        with open("tmp/" + region + "_profile.hmm", 'w') as profile_path:

            profile_path.write(query.profile.decode('utf-8'))



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


class GetUniqueSpecies(BaseSQLAFilter):

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
        print (species_list)
        print (id_list)
        return query.filter(self.get_column(alias).in_(id_list))

    def operation(self):
        return 'Yes'

class GetUniqueSpeciesSequenceRecord(BaseSQLAFilter):

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
                print(species_list)
                print(id_list)
                return query.filter(self.get_column(alias).in_(id_list))

            def operation(self):
                return 'Yes'


class SequenceRecordsView(ModelView):
    create_modal = True
    edit_modal = True
    can_create = False
    can_view_details = True

    # @expose('/')
    # def index(self):
    #     # Here are the contents of your "contact" route function
    #     return self.render('sequence_records.html', model=models.SequenceRecords, session=servers.db.session,)
    #


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

    # column_filters = (
    #     'sequence',
    #     GetUniqueSpeciesSequenceRecord(
    #         models.SequenceRecords.species, 'Get unique species')
    # )



    @action('set_A1_reference', 'Set this sequence as the A1 reference')
    def action_set_A1_reference(self, ids):

        setSequenceAsReference(ids, "a1")

    @action('set_A2_reference', 'Set this sequence as the A2 reference')
    def action_set_A2_reference(self, ids):

        setSequenceAsReference(ids, "a2")

    @action('set_pore_reference', 'Set this sequence as the pore reference')
    def action_set_pore_reference(self, ids):

        setSequenceAsReference(ids, "pore")

    @action('set_chitinase_reference', 'Set this sequence as the chitinase reference')
    def action_set_chitinase_reference(self, ids):

        setSequenceAsReference(ids, "chitinase")


    @action('set_region1_reference', 'Set this sequence as the region 1 reference')
    def action_set_region1_reference(self, ids):

        setSequenceAsReference(ids, "region1")

    @action('set_region2_reference', 'Set this sequence as the region 2 reference')
    def action_set_region2_reference(self, ids):

        setSequenceAsReference(ids, "region2")

    @action('set_region3_reference', 'Set this sequence as the region 3 reference')
    def action_set_region3_reference(self, ids):

        setSequenceAsReference(ids, "region3")

    @action('set_region4_reference', 'Set this sequence as the region 4 reference')
    def action_set_region4_reference(self, ids):

        setSequenceAsReference(ids, "region4")


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
    column_list = ('name', 'species', 'strain', 'description', 'a1_ref', 'a2_ref', 'pore_ref', 'sequence', 'a1', 'a1_length', 'a1_loc', 'a2',
                   'a2_length', 'a2_loc', 'overlap', 'distance', 'pore', 'pore_length', 'pore_loc', 'pore_within_a2',
                   'chitinase', 'chitinase_length', 'chitinase_loc', 'chitinase_distance_from_a2', 'region1_ref', 'region2_ref',
                   'region3_ref', 'region4_ref', 'region1', 'region1_length', 'region1_loc', 'region2',
                   'region2_length', 'region2_loc', 'region3', 'region3_length', 'region3_loc', 'region4',
                   'region4_length', 'region4_loc')
    column_searchable_list = ['name', 'species', 'a1', 'a2', 'overlap']
    create_modal = True
    edit_modal = True
    can_create = False
    can_view_details = True
    page_size = 20


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

    def _chitinase_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.chitinase:
            return model.chitinase[:15] + "..."
        else:
            return model.chitinase

    def _pore_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.pore:
            return model.pore[:15] + "..."
        else:
            return model.pore

    def _region1description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.region1:
            return model.region1[:15] + "..."
        else:
            return model.region1
    def _region2description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.region2:
            return model.region2[:15] + "..."
        else:
            return model.region2

    def _region3description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.region3:
            return model.region3[:15] + "..."
        else:
            return model.region3
    def _region4description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.region1:
            return model.region4[:15] + "..."
        else:
            return model.region4

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
        'pore': _pore_formatter,
        'chitinase': _chitinase_formatter,
        'region1': _region1description_formatter,
        'region2': _region2description_formatter,
        'region3': _region3description_formatter,
        'region4': _region4description_formatter,

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
                      GetUniqueSpecies(
                          models.GenomeRecords.name, 'Get unique species'),
                      filters.IntGreaterFilter(models.GenomeRecords.distance, 'Distance greater than'),
                      filters.IntSmallerFilter(models.GenomeRecords.distance, 'Distance smaller than'),
                      filters.IntGreaterFilter(models.GenomeRecords.chitinase_distance_from_a2, 'Chitinase to A2 distance greater than'),
                      filters.IntSmallerFilter(models.GenomeRecords.chitinase_distance_from_a2, 'Chitinase to A2 distance smaller than')

                      )

    def after_model_delete(self, model):
        """
        Function to delete the record from the BioSQL database as well
        :param model: Record to delete
        :return:
        """
        del servers.bio_db[model.uid]
        servers.bio_server.commit()

    @action('item5_get_closest_chitinase', 'Get closest chitinase to A2 region')
    def item5_get_closest_chitinase(self, ids):
        try:
            for id in ids:

                query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(id))
                for record in query.all():
                    if record.a2 and record.a2_loc:
                        a2_end_position = record.a2_loc.split(":")[1]
                        checkForRegion(id, 'chitinase', a2_end_position)



                    else:
                        flash("That record doesn't have an A2 region associated with it")

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to get closest chitinase %(error)s', error=str(ex)), 'error')

    @action('item4_check_if_pore_within_A2', 'Check if the pore is within the A2 region')

    def item4_check_if_pore_within_A2(self, ids):
        try:
            query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
            for record in query.all():
                if record.a2 == "" or record.pore == "":
                    pass
                else:
                    contains = phyloisland.check_if_contains(record.a2_loc, record.pore_loc)
                    if contains:
                        setattr(record, "pore_within_a2", "True")
                    else:
                        setattr(record, "pore_within_a2", "False")

                    servers.db.session.add(record)
                    servers.db.session.commit()

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to get overlap %(error)s', error=str(ex)), 'error')

    @action('item5_get_distance_between_A2_and_chitinase', 'Get distance between A2 and chitinase')
    def item5_get_distance_between_A2_and_chitinase(self, ids):
        try:
            query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
            for record in query.all():
                if record.a2 == "" or record.chitinase == "":
                    pass
                else:
                    distance = str(phyloisland.getDistance(record.a2_loc, record.chitinase_loc))
                    record.chitinase_distance_from_a2 = distance

                    servers.db.session.add(record)
                    servers.db.session.commit()

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to get overlap %(error)s', error=str(ex)), 'error')


    @action('item5_get_overlap', 'Get overlap')
    def item5_get_overlap(self, ids):
        try:
            query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
            for record in query.all():
                if record.a1 == "" or record.a2 == "":
                    pass
                else:
                    distance = str(phyloisland.getDistance(record.a1_loc, record.a2_loc))
                    print(distance)
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

        setRegionAsReference(ids, 'a1')


    @action('item7_set_A2_reference', 'Set this A2 region as the reference region')
    def item7_set_a2_reference(self, ids):

        setRegionAsReference(ids, 'a2')

    @action('item7_set_pore_reference', 'Set this pore region as the reference region')
    def item7_set_pore_reference(self, ids):

        setRegionAsReference(ids, 'pore')


    @action('item7_set_chitinase_reference', 'Set this chitinase region as the reference region')
    def item7_set_chitinase_reference(self, ids):

        setRegionAsReference(ids, 'chitinase')

    @action('item7_set_region1_reference', 'Set this region1 region as the reference region')
    def item7_set_region1_reference(self, ids):

        setRegionAsReference(ids, 'region1')

    @action('item7_set_region2_reference', 'Set this region2 region as the reference region')
    def item7_set_region2_reference(self, ids):

        setRegionAsReference(ids, 'region2')

    @action('item7_set_region3_reference', 'Set this region3 region as the reference region')
    def item7_set_region3_reference(self, ids):

        setRegionAsReference(ids, 'region3')

    @action('item7_set_region4_reference', 'Set this region4 region as the reference region')
    def item7_set_region4_reference(self, ids):

        setRegionAsReference(ids, 'region4')

    @action('item8_a1_profile_align_a1', 'Check for A1 region with a profile')
    def item_a1_profile_align_a1(self, ids):
        try:

            checkWithProfile(ids, 'a1')

        # Does this exception raise correctly?
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for A1 region -  %(error)s', error=str(ex)), 'error')


    @action('item8_a2_profile_align_a2', 'Check for A2 region with a profile')
    def item_a2_profile_align_a2(self, ids):
        try:
            checkWithProfile(ids, 'a2')
            # Does this exception raise correctly?
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for A2 region -  %(error)s', error=str(ex)), 'error')

    @action('item8_a2_profile_align_pore', 'Check for pore with a profile')
    def item_a2_profile_align_pore(self, ids):
        try:
            checkWithProfile(ids, 'pore')
            # Does this exception raise correctly?
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for pore -  %(error)s', error=str(ex)), 'error')

    @action('item8_a2_profile_align_chitinase', 'Check for chitinase with a profile')
    def item_a2_profile_align_chitinase(self, ids):
        try:
            checkWithProfile(ids, 'chitinase')
            # Does this exception raise correctly?
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for chitinase -  %(error)s', error=str(ex)), 'error')

    @action('item8_a2_profile_align_region1', 'Check for region 1 with a profile')
    def item_a2_profile_align_region1(self, ids):
        try:
            checkWithProfile(ids, 'region1')
            # Does this exception raise correctly?
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for region 1 -  %(error)s', error=str(ex)), 'error')

    @action('item_a2_profile_align_region2', 'Check for region 2 with a profile')
    def item_a2_profile_align_region2(self, ids):
        try:
            checkWithProfile(ids, 'region2')
            # Does this exception raise correctly?
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for region 2 -  %(error)s', error=str(ex)), 'error')

    @action('item_a2_profile_align_region3', 'Check for region 3 with a profile')
    def item_a2_profile_align_region3(self, ids):
        try:
            checkWithProfile(ids, 'region3')
            # Does this exception raise correctly?
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for region 3 -  %(error)s', error=str(ex)), 'error')

    @action('item_a2_profile_align_region4', 'Check for region 4 with a profile')
    def item_a2_profile_align_region4(self, ids):
        try:
            checkWithProfile(ids, 'region4')
            # Does this exception raise correctly?
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for region 4 -  %(error)s', error=str(ex)), 'error')

    @action('item3_check_a1', 'Check for A1 region')
    def item3_check_a1(self, ids):
        try:

            checkForRegion(ids, 'a1')


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for A1 region -  %(error)s', error=str(ex)), 'error')

    @action('item3_check_a2', 'Check for A2 region')
    def item3_check_a2(self, ids):
        try:

            checkForRegion(ids, 'a2')


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for A2 region -  %(error)s', error=str(ex)), 'error')

    @action('item3_check_pore', 'Check for pore')
    def item3_check_pore(self, ids):
        try:

            checkForRegion(ids, 'pore')


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for pore -  %(error)s', error=str(ex)), 'error')


    @action('item3_check_chitinase', 'Check for chitinase')
    def item3_check_chitinase(self, ids):
        try:

            checkForRegion(ids, 'chitinase')


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for chitinase -  %(error)s', error=str(ex)), 'error')

    @action('item3_check_region1', 'Check for region1')
    def item3_check_region1(self, ids):
        try:

            checkForRegion(ids, 'region1')


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for region 1 -  %(error)s', error=str(ex)), 'error')

    @action('item3_check_region2', 'Check for region2')
    def item3_check_region2(self, ids):
        try:

            checkForRegion(ids, 'region2')

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for region 2 -  %(error)s', error=str(ex)), 'error')

    @action('item3_check_region3', 'Check for region3')
    def item3_check_region3(self, ids):
        try:

            checkForRegion(ids, 'region3')


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for region 3 -  %(error)s', error=str(ex)), 'error')

    @action('item3_check_region4', 'Check for region4')
    def item3_check_region4(self, ids):
        try:

            checkForRegion(ids, 'region4')


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Something went wrong when checking for region 4 -  %(error)s', error=str(ex)), 'error')

    @action('item1_delete_a1', 'Delete A1 region')
    def item1_delete_a1(self, ids):
        try:

            checkForFeature.deleteFeature(ids=ids, string=["a1", "a1_loc"], int=["a1_length"], bool=["a1_ref"])

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to delete region. %(error)s', error=str(ex)), 'error')


    @action('item2_delete_a2', 'Delete A2 region')
    def item2_delete_a2(self, ids):
        try:
            checkForFeature.deleteFeature(ids=ids, string=["a2", "a2_loc", "pore_within_a2"], int=["a2_length"], bool=["a2_ref"])


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to delete region. %(error)s', error=str(ex)), 'error')

    @action('item2_delete_pore', 'Delete pore')
    def item2_delete_pore(self, ids):
        try:

            checkForFeature.deleteFeature(ids=ids, string=["pore", "pore_loc", "pore_within_a2"], int=["pore_length"], bool=["pore_ref"])


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to delete pore. %(error)s', error=str(ex)), 'error')

    @action('item2_delete_chitinase', 'Delete chitinase')
    def item2_delete_chitinase(self, ids):
        try:

            checkForFeature.deleteFeature(ids=ids, string=["chitinase", "chitinase_loc"], int=["chitinase_length", "chitinase_distance_from_a2"], bool=["chitinase_ref"])


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to delete chitinase. %(error)s', error=str(ex)), 'error')

    @action('item2_delete_region1', 'Delete region 1')
    def item2_delete_region1(self, ids):
        try:

            checkForFeature.deleteFeature(ids=ids, string=["region1", "region1_loc"], int=["region1_length"], bool=["region1_ref"])


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to delete region. %(error)s', error=str(ex)), 'error')
    @action('item2_delete_region2', 'Delete region 2')
    def item2_delete_region2(self, ids):
        try:

            checkForFeature.deleteFeature(ids=ids, string=["region2", "region2_loc"], int=["region2_length"], bool=["region2_ref"])

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to delete region. %(error)s', error=str(ex)), 'error')

    @action('item2_delete_region3', 'Delete region 3')
    def item2_delete_region3(self, ids):
        try:

            checkForFeature.deleteFeature(ids=ids, string=["region3", "region3_loc"], int=["region3_length"], bool=["region3_ref"])

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to delete region. %(error)s', error=str(ex)), 'error')

    @action('item2_delete_region4', 'Delete region 4')
    def item2_delete_region4(self, ids):
        try:

            checkForFeature.deleteFeature(ids=ids, string=["region4", "region4_loc"], int=["region4_length"], bool=["region4_ref"])

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise

            flash(gettext('Failed to delete region. %(error)s', error=str(ex)), 'error')

    @action('itema9_generate_profile_a1', 'Generate profile based on A1')
    def itema9_generate_profile_a1(self, ids):
        try:
            createProfileFromRegion(ids, 'a1')

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to generate profile based on A1. %(error)s', error=str(ex)), 'error')

    @action('itema9_generate_profile_a2', 'Generate profile based on A2')
    def itema9_generate_profile_a2(self, ids):
        try:
            createProfileFromRegion(ids, 'a2')
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to generate profile based on A2. %(error)s', error=str(ex)), 'error')

    @action('item9_generate_profile_pore', 'Generate profile based on pore')
    def itema9_generate_profile_pore(self, ids):
        try:
            createProfileFromRegion(ids, 'pore')
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to generate profile based on pore. %(error)s', error=str(ex)), 'error')

    @action('item9_generate_profile_chitinase', 'Generate profile based on chitinase')
    def itema9_generate_profile_chitinase(self, ids):
        try:
            createProfileFromRegion(ids, 'chitinase')
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to generate profile based on chitinase. %(error)s', error=str(ex)), 'error')

    @action('item9_generate_profile_region1', 'Generate profile based on region 1')
    def itema9_generate_profile_region1(self, ids):
        try:
            createProfileFromRegion(ids, 'region1')
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to generate profile based on region 1. %(error)s', error=str(ex)), 'error')

    @action('item9_generate_profile_region2', 'Generate profile based on region 2')
    def itema9_generate_profile_region2(self, ids):
        try:
            createProfileFromRegion(ids, 'region2')
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to generate profile based on region 2. %(error)s', error=str(ex)), 'error')

    @action('item9_generate_profile_region3', 'Generate profile based on region 3')
    def itema9_generate_profile_region3(self, ids):
        try:
            createProfileFromRegion(ids, 'region3')
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to generate profile based on region 3. %(error)s', error=str(ex)), 'error')

    @action('item9_generate_profile_region4', 'Generate profile based on region 4')
    def itema9_generate_profile_region4(self, ids):
        try:
            createProfileFromRegion(ids, 'region4')
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to generate profile based on region 4. %(error)s', error=str(ex)), 'error')

    @action('item9a_download_A2_sequences', 'Download the A2 sequences as a FASTA file')
    def item9a_download_A2_sequences(self, ids):
        createFASTAFromRegion(ids, "a2")



class ProfileView(ModelView):
    """
    View of the Profile database for storing HMM Profiles """
    column_list = ('name', 'a1_profile_ref', 'a2_profile_ref', 'pore_profile_ref', 'chitinase_profile_ref', 'region1_profile_ref', 'region2_profile_ref',
                   'region3_profile_ref', 'region4_profile_ref', 'download')
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

        setProfileAsReference(ids, "a1")


    @action('item2_set_A2_reference', 'Set this profile as the A2 reference profile')
    def item1_set_a2_reference(self, ids):
        setProfileAsReference(ids, "a2")

    @action('item2_set_pore_reference', 'Set this profile as the pore reference profile')
    def item1_set_pore_reference(self, ids):
        setProfileAsReference(ids, "pore")

    @action('item2_set_chitinase_reference', 'Set this profile as the chitinase reference profile')
    def item1_set_chitinase_reference(self, ids):
        setProfileAsReference(ids, "chitinase")

    @action('item2_set_region1_reference', 'Set this profile as the region 1 reference profile')
    def item1_set_region1_reference(self, ids):
        setProfileAsReference(ids, "region1")

    @action('item2_set_region2_reference', 'Set this profile as the region 2 reference profile')
    def item1_set_region2_reference(self, ids):
        setProfileAsReference(ids, "region2")

    @action('item2_set_region3_reference', 'Set this profile as the region 3 reference profile')
    def item1_set_region3_reference(self, ids):
        setProfileAsReference(ids, "region3")

    @action('item2_set_region4_reference', 'Set this profile as the region 4 reference profile')
    def item1_set_region4_reference(self, ids):
        setProfileAsReference(ids, "region4")


class UploadView(BaseView):
    """
    View for uploading files
    """

    @expose("/", methods =('GET', 'POST'))
    def upload(self):
        form = models.UploadForm()

        if request.method == 'POST' and form.validate():

            # Get the information from the upload form
            filename = servers.allfiles.save(request.files['file'])
            type = form.type.data
            add_seq = form.add_sequence.data
            add_genome = form.add_genome.data
            search_shotgun = form.search_shotgun.data

            # # Create the initial seqRecords
            # phyloisland.seqDict = {}
            # phyloisland.unmappable = []

            if type == "protein" or type == "nucleotide":
                seq_records = SeqIO.to_dict(SeqIO.parse("static/uploads/" + filename, "fasta"))

                if not seq_records:
                    print("Couldn't find any sequences in the uploaded file.")
                    flash("Couldn't find any sequences in the uploaded file.")

                else:
                    if add_seq:
                        addSequence(seq_records)

                    if add_genome:
                        species_names = mapToGenome.getSpeciesNames(seq_records, type)
                        genome_ids = mapToGenome.getGenomeIDs(species_names)
                        genome_query_string = utilities.makeQueryString(genome_ids, link="+OR+")

                        if genome_query_string == "":
                            print("We didn't identify any genome records. Attempting to search for shotgun sequenced "
                                  "genomes \n")
                            genome_results = mapToGenome.get_shotgun_id_dict(species_names)
                        else:
                            genome_results = mapToGenome.getFullGenome(genome_ids)

                        if genome_results and genome_results != 'in_database':
                            addGenome(genome_results)
                        elif genome_results != 'in_database' and search_shotgun:
                            print ("All of the genome records we identifed were all N characters. Attempting to "
                                   "search for shotgun sequenced genomes \n")
                            shotgun_id_dict = mapToGenome.get_shotgun_id_dict(species_names)
                            genome_results = mapToGenome.getFullGenome(shotgun_id_dict)

                            if genome_results:
                                addGenome(genome_results)

                        elif genome_results != 'in_database' and not search_shotgun:
                            print(
                                "\nWe didn't identify any genome records for %s. And we are not attempting to search for shotgun sequenced genomes \n" % (
                                    name))

            elif type == "species":
                species_names = readLinesFromFile("static/uploads/" + filename)

                for name in species_names:
                    genome_ids = {}
                    genome_results = {}
                    genome_ids = mapToGenome.getGenomeIDs(name)
                    genome_results = mapToGenome.getFullGenome(genome_ids)
                    print ('here is the genome results')
                    print (genome_results)
                    if genome_results and genome_results != 'in_database':
                        addGenome(genome_results)
                    elif genome_results != 'in_database' and search_shotgun:
                        print("\nWe didn't identify any genome records for %s. Attempting to search for shotgun sequenced genomes \n" % (name))
                        shotgun_id_dict = mapToGenome.get_shotgun_id_dict(name)
                        genome_results = mapToGenome.get_shotgun_genome(shotgun_id_dict)

                        if genome_results:
                            addGenome(genome_results)
                    elif genome_results != 'in_database' and not  search_shotgun:
                        print(
                            "\nWe didn't identify any genome records for %s. And we are not attempting to search for shotgun sequenced genomes \n" % (
                                name))


            elif (type == "genome"):
                genome_names = readLinesFromFile("static/uploads/" + filename)

                for name in genome_names:
                    genome_ids = {}
                    genome_results = {}
                    genome_ids = readLinesFromFile("static/uploads/" + filename)
                    genome_results = mapToGenome.getFullGenome(genome_ids)
                    if genome_results and genome_results != 'in_database':
                        addGenome(genome_results)
                    elif genome_results != 'in_database' and search_shotgun:
                        print(
                            "\nWe didn't identify any genome records for %s. Attempting to search for shotgun sequenced genomes \n" % (
                            name))

                        shotgun_id_dict = mapToGenome.get_shotgun_id_dict(name)
                        genome_results = mapToGenome.get_shotgun_genome(shotgun_id_dict)

                        if genome_results:
                            addGenome(genome_results)
                    elif genome_results != 'in_database' and not  search_shotgun:
                        print(
                            "\nWe didn't identify any genome records for %s. And we are not attempting to search for shotgun sequenced genomes \n" % (
                                name))


        return self.render("upload_admin.html", form=form)



class MyHomeView(AdminIndexView):
    @expose(servers.base_route)
    def index(self):
        return self.render('admin/index.html')


def addGenome(genome_results):

    for record in genome_results:

        current = genome_results[record]
        name = current.id

        species = " ".join(current.annotations.get('organism').split()[0:2])
        strain = current.annotations['source']
        sequence = str(current.seq)

        description = current.description
        a1 = current.annotations["A1"] if "A1" in current.annotations.keys() else ""
        a1_length = None
        a1_loc = current.annotations["A1_location"] if "A1_location" in current.annotations.keys() else ""
        a2 = current.annotations["A2"] if "A2" in current.annotations.keys() else ""
        a2_length = None
        a2_loc = current.annotations["A2_location"] if "A2_location" in current.annotations.keys() else ""
        overlap = current.annotations["Overlap"] if "Overlap" in current.annotations.keys() else ""
        distance = ""

        entry = models.GenomeRecords(name, species, strain, description, a1, a1_length, a1_loc, a2, a2_length, a2_loc,
                                     overlap, distance, sequence, 0, 0)
        check = models.GenomeRecords.query.filter_by(name=name).first()

        # Check to see if the genome record already exists
        if check:
            print("The genome record - %s from species - %s already exists in the database" % (name, species))

            continue
        # # if region in current.annotations.keys():
        #     #     continue
        #     # else:
        #     #     setattr(check, region.lower(), current.annotations[region] if region in current.annotations.keys() else "")
        #     #     setattr(check, region.lower() + "_loc", current.annotations[region + "_location"] if region + "_location" in current.annotations.keys() else "")
        #     #     servers.db.session.add(check)

        else:
            print("Adding the genome record - %s from species - %s to the genome database" % (name, species))

            seq_list = []
            servers.db.session.add(entry)
            seq_list.append(current)
            servers.bio_db.load(seq_list)
            servers.bio_server.commit()
            servers.db.session.commit()



def addSequence(seq_records):
    for record in seq_records.values():
        seq_name = record.id
        seq_description = record.description.split(">")[0]
        seq_species = seq_description.split("[")[1].split("]")[0]
        seq_sequence = str(record.seq)

        seq_entry = models.SequenceRecords(seq_name, seq_species, seq_description, seq_sequence, 0, 0)

        # Check if the sequence record already exists
        seq_check = models.SequenceRecords.query.filter_by(name=seq_name).first()
        if seq_check == None:
            print('Adding sequence with ID - %s from species - %s to the sequence database' % (
            seq_name, seq_species) + "\n")
            servers.db.session.add(seq_entry)
            servers.db.session.commit()
        else:
            print('Sequence with ID - %s from species - %s already exists in the sequence database' % (
            seq_name, seq_species) + "\n")

def setRegionAsReference(ids, region):
    if len(ids) > 1:
        flash('Only select a single record')
    else:
        query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
        for record in query.all():
            if (eval('record.' + region)):

                # Check for a previous reference sequence
                old_genome_reference = eval("models.GenomeRecords.query.filter_by(" + region + "_ref=1).first()")
                old_seq_reference = eval("models.SequenceRecords.query.filter_by(" + region + "_ref=1).first()")

                if (old_genome_reference):
                    # Remove the previous reference sequence
                    setattr(old_genome_reference, region + "_ref", 0)
                    servers.db.session.add(old_genome_reference)

                if (old_seq_reference):
                    # Remove the previous reference sequence
                    setattr(old_seq_reference, region + "_ref", 0)
                    servers.db.session.add(old_seq_reference)

                # Set the new reference sequence
                setattr(record, region + "_ref", 1)

                # Commit the changed record
                servers.db.session.add(record)
                servers.db.session.commit()
                flash("The %s region from %s has been set as the reference %s sequence" % (region, record.name, region))
            else:
                flash("That record doesn't have %s annotated onto it" % (region))


def setSequenceAsReference(ids, region):
    if len(ids) > 1:
        flash('Only select a single record')
    else:
        query = models.SequenceRecords.query.filter(models.SequenceRecords.uid.in_(ids))
        for record in query.all():

            # Check for a previous reference sequence
            old_genome_reference = eval("models.GenomeRecords.query.filter_by(" + region + "_ref=1).first()")
            old_seq_reference = eval("models.SequenceRecords.query.filter_by(" + region + "_ref=1).first()")

            if old_genome_reference:
                # Remove the previous reference sequence
                setattr(old_genome_reference, region + "_ref", 0)
                servers.db.session.add(old_genome_reference)

            if old_seq_reference:
                # Remove the previous reference sequence
                setattr(old_seq_reference, region + "_ref", 0)
                servers.db.session.add(old_seq_reference)

            # Set the new reference sequence
            setattr(record, region + "_ref", 1)

            # Commit the changed record
            servers.db.session.add(record)
            servers.db.session.commit()

            # global a1_reference
            region_name = region + "_reference"
            # eval("region_name") = SeqRecord(Seq(record.sequence, generic_protein), id=str(record.name) + '_reference')
            flash("The %s region from %s has been set as the reference %s sequence" % ( region, record.name, region))

def setProfileAsReference(ids, region):
    if len(ids) > 1:
        flash('Only select a single record')
    else:
        query = models.Profile.query.filter(models.Profile.uid.in_(ids))
        for record in query.all():

            # Check for a previous reference profile
            old_profile_reference = eval("models.Profile.query.filter_by(" + region + "_profile_ref=1).first()")

            if (old_profile_reference):
                # Remove the previous reference profile
                setattr(old_profile_reference, region + "_profile_ref", 0)
                servers.db.session.add(old_profile_reference)

            # Set the new reference profile
            setattr(record, region + "_profile_ref", 1)

            # Commit the changed record
            servers.db.session.add(record)
            servers.db.session.commit()

            # Write the new profile to the tmp folder ready to be used

            with open("tmp/" + region + "_profile.hmm", 'w') as profile_path:
                profile_path.write(record.profile.decode('utf-8'))

            flash("The profile named %s has been set as the reference profile for %s" % (record.name, region))

def checkForRegion(ids, region, closest_to = -1):
    genome_reference = eval("models.GenomeRecords.query.filter_by(" + region + "_ref=1).first()")
    sequence_reference = eval("models.SequenceRecords.query.filter_by(" + region + "_ref=1).first()")



    if genome_reference:
        reference = eval('SeqRecord(Seq(genome_reference.' + region + ', generic_protein), id=str(genome_reference.name) + "_" + "' + region + '")')
        checkForFeature.getFeatureLocation(ids, reference, region, region + "_loc", region + "_length", closest_to)


    elif sequence_reference:
        reference = eval('SeqRecord(Seq(sequence_reference.sequence, generic_protein), id=str(sequence_reference.name) + "_" + "' + region + '")')
        checkForFeature.getFeatureLocation(ids, reference, region, region + "_loc", region + "_length", closest_to)


    else:
        flash ("Please set a reference sequence for %s first" % (region))



def checkWithProfile(ids, region):

    # Check if a reference profile for A1 exists
    profile_reference = eval("models.Profile.query.filter_by(" + region + "_profile_ref=1).first()")

    if (profile_reference):
        print("Using the %s profile named %s to check for %s regions" % (region, profile_reference.name, region))
        eval('checkForFeature.get_feature_location_with_profile(ids, "tmp/' + region + '_profile.hmm", "' + region + '", "' + region + ' _loc")')
    else:
        flash("Please set a profile as the %s reference profile first" % (region), "error")



def createProfile(align_list):
    SeqIO.write(align_list, "align.fasta", "fasta")
    muscle_cline = MuscleCommandline(input="align.fasta")
    # result = subprocess.run(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True) )
    child = subprocess.Popen(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True, shell=(sys.platform != "win32"))
    child.wait()

    alignment = AlignIO.read(child.stdout, "fasta")
    AlignIO.write(alignment, "align.aln", "fasta")
    hmm_path = "tmp/profile3.hmm"
    outfile = open(hmm_path, "w")
    result = subprocess.call(["hmmbuild", hmm_path, "align.aln"], stdout=subprocess.PIPE)


    while not os.path.exists(hmm_path):
        time.sleep(1)

    if os.path.isfile(hmm_path):

        file = open(hmm_path, 'rb')


        saveProfile(file)
        utilities.removeFile(hmm_path, "align.fasta", "align.aln")

def createProfileFromRegion(ids, region):
        query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
        align_list = []
        for record in query.all():
            if eval('record.' + region  + '== ""'):
                pass
            else:
                align_record = eval('SeqRecord(Seq(record.' + region + ', generic_protein), id=str(record.name) + "_" + "' + region +'")')
                align_list.append(align_record)

        if align_list:
            createProfile(align_list)
        else:
            flash("None of the selected genomes has an annotated %s so we couldn't build a profile" % (region))


def createFASTAFromRegion(ids, region):
    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
    align_list = []
    for record in query.all():
        if eval('record.' + region + '== ""'):
            pass
        else:
            align_record = eval(
                'SeqRecord(Seq(record.' + region + ', generic_protein), id=str(record.name), description="' + region + '_region")')
            align_list.append(align_record)

    if align_list:
        utilities.saveFASTA(align_list, "tmp/output.fasta")
    else:
        flash("None of the selected genomes has an annotated %s so we couldn't build a profile" % (region))



def saveProfile(profile):
    """
    Save a profile into the database
    :param profile: Profile to save
    """

    name = phyloisland.randstring(5)
    blobMix = models.BlobMixin("application/octet-stream", name, profile.read(), '566666')
    profileEntry = models.Profile(name)
    profileEntry.set_blobMix(blobMix)
    servers.db.session.add(profileEntry)
    servers.db.session.commit()

def readLinesFromFile(filepath):
    """
    Takes a file and reads each individual line into a set
    :param filepath: Path of the file
    :return: Set containing lines from the file
    """

    content = set()

    with open(filepath, 'r') as query_file:
        for line in query_file:
            if len(line) > 1:
                content.add(line.strip())
    return content

# Setup the main flask-admin application
admin = Admin(servers.application, index_view=AdminIndexView(name= servers.base_name, url="/" + servers.base_route), template_mode="bootstrap3")

# Add the views to the flask-admin application
admin.add_view(UploadView(name='Upload', endpoint='upload_admin'))
admin.add_view(SequenceRecordsView(model=models.SequenceRecords, session=servers.db.session, endpoint="sequence_records"))
admin.add_view(GenomeRecordsView(model=models.GenomeRecords, session=servers.db.session, endpoint="genome_view"))
admin.add_view(ProfileView(model=models.Profile, session=servers.db.session, name='Profiles'))

# # Create A1 and A2 reference profiles on the disk
# setA1ReferenceProfile('a1')
# setA2ReferenceProfile('a2')

# Create A1 and A2 reference profiles on the disk
for region in ["a1", "a2", "pore", "chitinase", "region1", "region2", "region3", "region4"]:
    setReferenceProfile(region)


if __name__ == '__main__':
    servers.application.run(debug=True, port=7777, use_reloader=False)

