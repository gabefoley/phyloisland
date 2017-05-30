# import os
# from flask_table import Table, Col, LinkCol
# from flask import Flask, Markup, request, url_for, render_template, jsonify, redirect, flash, send_from_directory
# from Table import SortableTable, Item
# from Bio import SeqIO
# from forms import UploadForm
# from werkzeug.utils import secure_filename
#
# import subprocess
# from datetime import datetime
# """
# A example for creating a Table that is sortable by its header
# """
# from os.path import join, dirname, realpath
#
# UPLOADS_PATH = join(dirname(realpath(__file__)), 'static/uploads/..')
#
# app = Flask(__name__)
#
# APP_ROOT = os.path.dirname(os.path.abspath(__file__))
# app.config.from_object(__name__)
# app.config['UPLOAD_FOLDER'] = UPLOADS_PATH
#
# app.secret_key = 'development-key'
#
#
#
# # @app.route('/')
# # def index():
# #     form = UploadForm()
# #     sort = request.args.get('sort', 'id')
# #     reverse = (request.args.get('direction', 'asc') == 'desc')
# #     table = SortableTable(Item.get_sorted_by(sort, reverse),
# #                           sort_by=sort,
# #                           sort_reverse=reverse)
# #     return render_template("index.html", table=table, form = form)
#
#
# @app.route('/', methods=['GET', 'POST'])
# def upload_file():
#     sort = request.args.get('sort', 'id')
#     reverse = (request.args.get('direction', 'asc') == 'desc')
#     table = SortableTable(Item.get_sorted_by(sort, reverse),
#                               sort_by=sort,
#                               sort_reverse=reverse)
#     form = UploadForm()
#
#     if request.method == 'POST':
#         # check if the post request has the file part
#         if 'file' not in request.files:
#             flash('No file part')
#             return redirect(request.url)
#         file = request.files['file']
#         # if user does not select file, browser also
#         # submit a empty part without filename
#         if file.filename == '':
#             flash('No selected file')
#             return redirect(request.url)
#         if file and allowed_file(file.filename):
#             filename = secure_filename(file.filename)
#             file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
#             return redirect(url_for('uploaded_file', filename=filename))
#             # return render_template("index.html")
#     return render_template("index.html", form = form)
#
# @app.route('/uploads/<filename>')
# def uploaded_file(filename):
#     return send_from_directory(app.config['UPLOAD_FOLDER'],
#                                filename)
# # @app.route('/', methods=['GET', 'POST'])
# # def upload():
# #     sort = request.args.get('sort', 'id')
# #     reverse = (request.args.get('direction', 'asc') == 'desc')
# #     table = SortableTable(Item.get_sorted_by(sort, reverse),
# #                               sort_by=sort,
# #                               sort_reverse=reverse)
# #     form = UploadForm()
# #
# #     if form.validate_on_submit():
# #         # filename1 = secure_filename(form.file1.data.filename)
# #         # form.file1.data.save('uploads/' + filename1)
# #         # filename2 = secure_filename(form.file2.data.filename)
# #         # form.file2.data.save('uploads/' + filename2)
# #         # return redirect(url_for('upload'))
# #
# #         return render_template('index.html', form = form, table=table)
# #
# #     return render_template('index.html', form = form)
#
#
#     # if request.method == 'POST':
#     #     file = request.files['file']
#     #     if file and allowed_file(file.filename):
#     #         now = datetime.now()
#     #         filename = os.path.join(app.config['UPLOAD_FOLDER'], "%s.%s" % (now.strftime("%Y-%m-%d-%H-%M"), file.filename.rsplit('.', 1)[1]))
#     #         file.save(filename)
#     #         sq = filename
#     #         position = -1
#     #         for seq_record in SeqIO.parse(sq, "fasta"):
#     #             x = request.form['search_item']
#     #             seq = Seq(x)
#     #             seq_rev = seq.reverse_complement()
#     #             data = seq_record.seq
#     #             position = data.find(x)
#     #             position_rev = data.find(seq_rev)
#     #         return render_template('result.html', position=position, position_rev=position_rev, form = form)
#     #     else:
#     #         return render_template("no_fasta.html")
#
# def allowed_file(filename):
#     return filename
#
# @app.route('/item/<int:id>')
# def flask_link(id):
#     element = Item.get_element_by_id(id)
#     return '<h1>{}</h1><p>{}</p><hr><small>id: {}</small>'.format(
#         element.name, element.description, element.id)
#
#
#
# if __name__ == '__main__':
#     app.run(debug=True)

from flask_table import Table, Col, LinkCol
from flask import Flask, Markup, request, url_for
import Table

"""
A example for creating a Table that is sortable by its header
"""

app = Flask(__name__)


class SortableTable(Table):
    id = Col('ID')
    name = Col('Name')
    description = Col('Description')
    link = LinkCol(
        'Link', 'flask_link', url_kwargs=dict(id='id'), allow_sort=False)
    allow_sort = True

    def sort_url(self, col_key, reverse=False):
        if reverse:
            direction = 'desc'
        else:
            direction = 'asc'
        return url_for('index', sort=col_key, direction=direction)


@app.route('/')
def index():
    sort = request.args.get('sort', 'id')
    reverse = (request.args.get('direction', 'asc') == 'desc')
    table = SortableTable(Item.get_sorted_by(sort, reverse),
                          sort_by=sort,
                          sort_reverse=reverse)
    return table.__html__()


@app.route('/item/<int:id>')
def flask_link(id):
    element = Item.get_element_by_id(id)
    return '<h1>{}</h1><p>{}</p><hr><small>id: {}</small>'.format(
        element.name, element.description, element.id)


class Item(object):
    """ a little fake database """
    def __init__(self, id, name, description):
        self.id = id
        self.name = name
        self.description = description

    @classmethod
    def get_elements(cls):
        return [
            Item(1, 'Z', 'zzzzz'),
            Item(2, 'K', 'aaaaa'),
            Item(3, 'B', 'bbbbb')]

    @classmethod
    def get_sorted_by(cls, sort, reverse=False):
        return sorted(
            cls.get_elements(),
            key=lambda x: getattr(x, sort),
            reverse=reverse)

    @classmethod
    def get_element_by_id(cls, id):
        return [i for i in cls.get_elements() if i.id == id][0]

if __name__ == '__main__':
    app.run(debug=True)