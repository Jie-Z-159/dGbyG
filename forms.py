from flask_wtf import FlaskForm
from wtforms import Form, StringField, SelectField, SubmitField
from wtforms.validators import DataRequired

class ReactionForm(FlaskForm):
    equation = StringField (
        'Equation',validators=[DataRequired()],
        render_kw={"placeholder": "Example: 2 [H][H] + O=O = 2 O"}
    )
    identifier_type = SelectField(
        'Identifier type',
        choices=[('SMILES', 'SMILES'),
                 ('InChI', 'InChI'),
                 ('KEGG ID', 'KEGG ID'),
                 ('MetaNetX ID', 'MetaNetX ID'),
                 ('HMDB ID', 'HMDB ID')],
        validators=[DataRequired()]
    )    
    cell_compartment = SelectField(
        'Choose a cell compartment',
        choices=[('d', 'Default'),
                 ('c', 'Cytosol'),
                 ('e', 'Extracellular'),
                 ('n', 'Nucleus'),
                 ('r', 'Endoplasmic Reticulum'),
                 ('g', 'Golgi Apparatus'),
                 ('l', 'Lysosome'),
                 ('m', 'Mitochondria'),
                 ('i', 'Inner Mitochondria'),
                 ('x', 'Peroxisome')],
        validators=[DataRequired()]
    )
    submit = SubmitField('prediction')