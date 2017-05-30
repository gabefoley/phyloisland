from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, SubmitField, TextAreaField, FileField
from wtforms.validators import DataRequired, Email, Length, ValidationError

class UploadForm(FlaskForm):
    name = StringField('What ID should we give this feature?', validators=[DataRequired("Not completed")])
    file = FileField('Upload the FASTA file', validators=[DataRequired("Not completed")] )

    upload_submit = SubmitField("Upload files")

