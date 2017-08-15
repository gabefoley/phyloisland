from BioSQL import BioSeqDatabase
from flask_uploads import UploadSet, configure_uploads, ALL
from flask_sqlalchemy import SQLAlchemy
from flask import Flask

bio_server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="fishtank")
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