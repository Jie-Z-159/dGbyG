from flask_wtf import FlaskForm
from wtforms import Form, StringField, SelectField, SubmitField,RadioField
from wtforms.validators import DataRequired

class ReactionForm(FlaskForm):
    equation = StringField (
        'Equation',validators=[DataRequired(message="Please enter the chemical reaction formula")],
        render_kw={"placeholder": "Example: 2 [H][H] + O=O = 2 O"}
    )
    identifier_type = SelectField(
        'Identifier type',
        choices=[
            ('smiles', 'SMILES'),
            ('inchi', 'InChI'),
            ('kegg', 'KEGG ID'),
            ('metanetx', 'MetaNetX ID'),
            ('hmdb', 'HMDB ID')
        ],
        validators=[DataRequired()]
    )
    reaction_condition = RadioField(
        'Choose reaction condition',
        choices=[('d', 'Default'),
                 ('c', 'Cytosol'),
                 ('e', 'Extracellular'),
                 ('n', 'Nucleus'),
                 ('r', 'Endoplasmic Reticulum'),
                 ('g', 'Golgi Apparatus'),
                 ('l', 'Lysosome'),
                 ('m', 'Mitochondria'),
                 ('i', 'Inner Mitochondria'),
                 ('x', 'Peroxisome'),
                 ('custom','other condition')],
        validators=[DataRequired()]
    )
    submit = SubmitField('prediction')