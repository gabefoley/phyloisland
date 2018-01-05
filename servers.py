from BioSQL import BioSeqDatabase
from flask_uploads import UploadSet, configure_uploads, ALL
from flask_sqlalchemy import SQLAlchemy
from flask import Flask

# Names to use for our parameters
base_name = "Experimental"
base_route = "phyloisland_experimental"
bio_server_name = "phyloisland"
bio_db_name = "newsmall"


# Setup the database
bio_server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db=bio_server_name)
bio_db = bio_server[bio_db_name]

# Setup the flask application
application = Flask(__name__)
application.config['SQLALCHEMY_DATABASE_URI'] = 'mysql://pi:@localhost/' + bio_server_name
application.config['SECRET_KEY'] = 'developmentkey'
application.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = True

allfiles = UploadSet('all', ALL)
application.config['UPLOADS_ALL_DEST'] = 'static/uploads'
application.config['UPLOADED_ALL_DEST'] = 'static/uploads'
configure_uploads(application, allfiles)

db = SQLAlchemy(application)