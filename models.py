from flask_sqlalchemy import SQLAlchemy
from werkzeug.security import generate_password_hash, check_password_hash
from admin import db

db = SQLAlchemy()

# class User(db.Model):
#     __tablename__ = 'users'
#     uid = db.Column(db.Integer, primary_key=True)
#     firstname = db.Column(db.String(100))
#     lastname = db.Column(db.String(100))
#     studentno = db.Column(db.Integer)
#     email = db.Column(db.String(120), unique=True)
#     pwdhash = db.Column(db.String(54))
#
#     def __init__(self, firstname="", lastname="", studentno="", email="", password=""):
#         self.firstname = firstname.title()
#         self.lastname = lastname.title()
#         self.studentno = studentno.title()
#         self.email = email.lower()
#         self.set_password(password)
#
#     def set_password(self, password):
#         self.pwdhash = generate_password_hash(password)
#
#     def check_password(self, password):
#         return check_password_hash(self.pwdhash, password)



class BioEntry(db.Model):
    __tablename__ = 'bioentry'
    bioentry_id = db.Column(db.Integer, primary_key=True)
    biodatabase_id = db.Column(db.Integer, unique = True)
    taxon_id = db.Column(db.Integer)
    name = db.Column(db.VARCHAR(40))
    accession = db.Column(db.VARCHAR(128), unique=True)
    identifier = db.Column(db.VARCHAR(40))
    division = db.Column(db.VARCHAR(6))
    description = db.Column(db.TEXT)
    version = db.Column(db.SMALLINT, unique=True)




class SeqRecord(db.Model):
    __tablename__ = 'seqrecord'
    uid = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(255))
    species = db.Column(db.String(255))
    strain = db.Column(db.String(255))
    description = db.Column(db.String(255))
    sequence =  db.Column(db.String(255))
