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
from flask.ext.admin.form import rules
from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, SubmitField, TextAreaField, FileField
from wtforms.validators import DataRequired, Email, Length, ValidationError
import subMenu
from flask_uploads import UploadSet, configure_uploads, ALL
from Table import SortableTable, Item
from flask_admin.actions import action
import gettext

subMenu.seqRecords = {}

# Create directory for file fields to use
file_path = op.join(op.dirname(__file__), 'files')
try:
    os.mkdir(file_path)
except OSError:
    pass

bio_server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="bioseqdb")
bio_db = bio_server["phylotest"]

app = Flask(__name__)
# app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql://pi:@localhost/bioseqdb'

app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql:///gabe'
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



class SeqRecord(db.Model):
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
    fullrecord = db.Column(db.Text)

    def __init__(self, name="", species="", strain="", description="", a1="", a1_loc="", a2="", a2_loc="", overlap="", distance="", sequence="", fullrecord=""):
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
        self.fullrecord = fullrecord

class SeqView(ModelView):
    column_sortable_list = ['name', 'division']

class SeqRecordView(ModelView):
    column_searchable_list = ['name', 'species', 'a1', 'a2']

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

    @action('getoverlap', 'Get Overlap')
    def action_approve(self, ids):
        print (ids)
        try:
            query = SeqRecord.query.filter(SeqRecord.uid.in_(ids))

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


class UploadForm(FlaskForm):
    name = StringField('What ID should we give this feature?', validators=[DataRequired("Not completed")])
    file = FileField('Upload the FASTA file', validators=[DataRequired("Not completed")] )

    upload_submit = SubmitField("Upload files")


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
        seqList = []

        if request.method == 'POST':
            # region = request.form['name']

            # Get the file and ID to name these sequences
            region = form.name.data
            filename = allfiles.save(request.files['file'])
            print (filename, region, "##")

            # Create the initial seqRecords
            subMenu.defaultValue("static/uploads/"  + filename, region)

            records = subMenu.seqRecords

            for record in records:
                current = records[record]
                # print (current)
                name = current.id
                species = current.annotations['source']
                strain = current.annotations['source']
                description = current.description
                a1 = current.annotations["A1"] if "A1" in current.annotations.keys() else "Not tested"
                a1_loc = current.annotations["A1_location"] if "A1_location" in current.annotations.keys() else "Not tested"
                a2 = current.annotations["A2"] if "A2" in current.annotations.keys() else "Not tested"
                a2_loc = current.annotations["A2_location"] if "A2_location" in current.annotations.keys() else "Not tested"
                overlap = current.annotations["Overlap"] if "Overlap" in current.annotations.keys() else "Not tested"
                distance = "Not tested"
                # sequence = str(current.seq)
                sequence = str(current.seq)
                fullrecord = current

                entry = SeqRecord(name, species, strain, description, a1, a1_loc, a2, a2_loc, overlap, distance, sequence)

                check = SeqRecord.query.filter_by(name=name).first()

                # check = db.session.query(db.exists().where(SeqRecord.name== name)).scalar()
                print ('Name we are checking on is ', name)

                print (check)

                if (check):
                    print ('got here baby')
                    print ("A1 loc ", check.a1_loc)
                    print ("A2 loc ", check.a2_loc)
                    check = subMenu.addToRecord(check, current, region)
                    print ("A1 loc ", check.a1_loc)
                    print ("A2 loc ", check.a2_loc)

                    db.session.add(check)
                else:
                    db.session.add(entry)
                db.session.commit()

            return self.render("upload_admin.html", form=form, records = records)
        return self.render("upload_admin.html", form=form)



admin = Admin(app, name="Phylo Island", template_mode="bootstrap3")
admin.add_view(UploadView(name='Upload', endpoint='upload_admin'))

admin.add_view(SeqRecordView(SeqRecord, db.session))
# admin.add_view(SeqView(BioEntry, db.session))
admin.add_view(UserView(User, db.session))



# admin.add_view(SeqRecordView(SeqRecord, db.session))



path = op.join(op.dirname(__file__), 'static')
admin.add_view(FileAdmin(path, '/static/', name='Static Files'))
# admin.add_view(FileView(File, db.session))
# admin.add_view(FileView(File, db.session))
# admin.add_view(UserView(User, db.session, name='User'))
# admin.add_view(AnalyticsView(name='Analytics', endpoint='analytics'))


@app.route('/')
def index():
    form = UploadForm()
    # return render_template("index.html", form=form)
    return '<a href="/admin/">Click me to get to Admin!</a> <a href="/upload">Click me to get to Upload!</a>'

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






app.run()

