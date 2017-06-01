import os
import subMenu
from flask import Flask, render_template, request, session, redirect, url_for, flash, send_from_directory
from werkzeug.utils import secure_filename
from forms import UploadForm
from flask_uploads import UploadSet, configure_uploads, ALL
from Table import SortableTable, Item
from flask_table import Table, Col, LinkCol
from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="bioseqdb")
db = server["phylotest"]
app = Flask(__name__)
app.secret_key = 'development-key'


# app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql:///gabe'


# db.init_app(app)
#
#
# id = 1
#
#
allfiles = UploadSet('all', ALL)
app.config['UPLOADS_ALL_DEST'] = 'static/uploads'

app.config['UPLOADED_ALL_DEST'] = 'static/uploads'
configure_uploads(app, allfiles)


# def allowed_file(filename):
#   return '.' in filename and filename.rsplit('.', 1)[1] in app.ognfig['ALLOWED_EXTENSIONS']

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
# app.secret_key = 'development-key'



@app.route("/", methods = ['GET', 'POST'])
def index():
    form = UploadForm()
    return render_template("index2.html", form=form)

@app.route("/upload", methods = ['GET', 'POST'])

def upload():
    form = UploadForm()
    sort = request.args.get('sort', 'id')
    reverse = (request.args.get('direction', 'asc') == 'desc')
    seqList = []

    if request.method == 'POST' and 'file' in request.files:
        count = 1
        item = Item(0, "Hard", "Coded")

        # Get the file and ID to name these sequences
        name = form.name.data
        filename = allfiles.save(request.files['file'])

        # Create the initial seqRecords
        subMenu.defaultValue("static/uploads/"  + filename, name)

        records = subMenu.seqDict

        for record in subMenu.seqDict:
            print (subMenu.seqDict[record].id)
            seqList.append(subMenu.seqDict[record])

            next_item = Item(count, subMenu.seqDict[record].id, subMenu.seqDict[record].annotations['organism'])
            item.add_elements(next_item)
            count += 1

            table = SortableTable(item.get_sorted_by(sort, reverse),
                                  sort_by=sort,
                                  sort_reverse=reverse)
        # session['item'] = item

        # Load into server
        # count = db.load(seqList)
        # print ('Loaded %d sequences' % count)
        # server.commit()

        # print (subMenu.seqRecords)
        return render_template("index.html", form=form, records = records)
    return render_template("index.html", form=form)

# @app.route('/item/<int:id>')
# def flask_link(id):
#     Item = session['item']
#     element = Item.get_element_by_id(id)
#     return '<h1>{}</h1><p>{}</p><hr><small>id: {}</small>'.format(
#       element.name, element.description, element.id)

record = db.lookup(accession = "CPYD01000004")


if __name__ == "__main__":
  app.run(debug=True, port = 5080)



  # @app.route('/')
  # def index():
  #     sort = request.args.get('sort', 'id')
  #     reverse = (request.args.get('direction', 'asc') == 'desc')
  #     table = SortableTable(Item.get_sorted_by(sort, reverse),
  #                           sort_by=sort,
  #                           sort_reverse=reverse)
  #     return table.__html__()











