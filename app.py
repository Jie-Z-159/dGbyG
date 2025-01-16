import os
import sys
import click
import logging
from flask import Flask, render_template, request,jsonify
from flask_sqlalchemy import SQLAlchemy
from dGbyG.api import Compound, Reaction
from dGbyG.utils.ChemFunc import parse_equation
from forms import ReactionForm

# 配置日志以调试
logging.basicConfig()
logging.getLogger('sqlalchemy.engine').setLevel(logging.INFO)

WIN = sys.platform.startswith('win')
if WIN:  # 如果是 Windows 系统，使用三个斜线
    prefix = 'sqlite:///'
else:  # 否则使用四个斜线
    prefix = 'sqlite:////'

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = prefix + os.path.join(app.root_path, 'data.db')
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.secret_key = 'your_secret_key_for_session'  # 请替换为你的密钥
db = SQLAlchemy(app)

class ReactionModel(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    identifier = db.Column(db.String(20))
    equation = db.Column(db.Text, nullable=False)
    reaction_condition = db.Column(db.String(20))
    dG_prime = db.Column(db.Float, nullable=True)
    dG_std_dev = db.Column(db.Float, nullable=True) 

    def __repr__(self):
        return '<ReactionModel %r>' % self.equation

# 创建数据库和表
@app.cli.command('initdb')
def initdb():
    db.create_all()
    click.echo('Initialized database.')

# 定义固定条件
default_T = 298.15
default_I = 0.25
default_pMg = 14.0


# 首页和预测路由
@app.route('/pasre_equation',methods = ['POST'])
def parse_equation_route():
    equation = request.form.get('equation')
    if not equation:
        return jsonify({'error': 'No equation provided'}), 400
    try:
        equation_dict = parse_equation(equation)
        compound_names = list(equation_dict.keys())
        return jsonify({'compounds': compound_names})
    except Exception as e:
        return jsonify({'error':str(e)}),500

@app.route('/', methods=['GET', 'POST'])
@app.route('/prediction', methods=['GET', 'POST'])
def prediction():
    form = ReactionForm()
    dG_prime = None
    dG_std_dev = None
    error_message = None
    conditions = {
        'd': {'pH': 7.00, 'e_potential': 0, 'T': default_T, 'I': default_I, 'pMg': default_pMg},
        'c': {'pH': 7.20, 'e_potential': 0, 'T': default_T, 'I': default_I, 'pMg': default_pMg},
        'e': {'pH': 7.40, 'e_potential': 30 * 1e-3, 'T': default_T, 'I': default_I, 'pMg': default_pMg},
        'n': {'pH': 7.20, 'e_potential': 0, 'T': default_T, 'I': default_I, 'pMg': default_pMg},
        'r': {'pH': 7.20, 'e_potential': 0, 'T': default_T, 'I': default_I, 'pMg': default_pMg},
        'g': {'pH': 6.35, 'e_potential': 0, 'T': default_T, 'I': default_I, 'pMg': default_pMg},
        'l': {'pH': 5.50, 'e_potential': 19 * 1e-3, 'T': default_T, 'I': default_I, 'pMg': default_pMg},
        'm': {'pH': 8.00, 'e_potential': -155 * 1e-3, 'T': default_T, 'I': default_I, 'pMg': default_pMg},
        'i': {'pH': 8.00, 'e_potential': -155 * 1e-3, 'T': default_T, 'I': default_I, 'pMg': default_pMg},
        'x': {'pH': 7.00, 'e_potential': 12 * 1e-3, 'T': default_T, 'I': default_I, 'pMg': default_pMg}
    }    

    if request.method == 'POST':
        if form.reaction_condition.data == 'custom':
            equation = form.equation.data
            identifier = form.identifier_type.data
            try:
                equation_dict = parse_equation(equation)
                compound_names = list(equation_dict.keys())

                custom_conditions = {}
                for index, compound in enumerate(compound_names):
                    pH = float(request.form.get(f'custom_conditions-{index}-pH', 7.0))
                    I = float(request.form.get(f'custom_conditions-{index}-I', 0.1))
                    pMg = float(request.form.get(f'custom_conditions-{index}-pMg', 7.0))
                    e_potential = float(request.form.get(f'custom_conditions-{index}-e_potential', 0.0))

                    custom_conditions[compound]= {
                        'pH': pH,
                        'I': I,
                        'pMg': pMg,
                        'e_potential': e_potential
                    }

                    # 创建化合物对象并设置自定义条件
                    compound_objs = []
                    stoichiometry = []
                    for compound_name in compound_names:
                        compound_obj = Compound(compound_name, input_type=identifier)
                        compound_obj.condition = custom_conditions.get(compound_name,{})
                        compound_objs.append(compound_obj)
                        stoichiometry.append(equation_dict[compound_name])
                    equation_dict_compounds = dict(zip(compound_objs, stoichiometry))
                    reaction = Reaction(equation_dict_compounds, cid_type='compound')
                    dG_prime, dG_std_dev = reaction.transformed_standard_dGr_prime

                    # 将数据保存到数据库
                    new_reaction = ReactionModel(
                        identifier=identifier,
                        equation=form.equation.data,
                        reaction_condition='custom',
                        dG_prime=dG_prime,
                        dG_std_dev=dG_std_dev
                    )
                    db.session.add(new_reaction)
                    db.session.commit()

            except Exception as e:
                error_message = str(e)
                logging.error(f"Error during prediction: {e}")

        else:
            if form.validate_on_submit():
                reaction_condition = form.reaction_condition.data
                identifier = form.identifier_type.data
                equation = form.equation.data

                try:
                    chosen_condition = conditions.get(reaction_condition)
                    reaction = Reaction(equation, cid_type=identifier)

                    for compound,coeff in reaction.reaction.items():
                        compound.condition['pH'] = chosen_condition['pH']
                        compound.condition['I'] = chosen_condition['I']
                        compound.condition['pMg'] = chosen_condition['pMg']
                        compound.condition['e_potential'] = chosen_condition['e_potential']
                    dG_prime, dG_std_dev = reaction.transformed_standard_dGr_prime

                    # 将数据保存到数据库
                    new_reaction = ReactionModel(
                        identifier=identifier,
                        equation=equation,
                        reaction_condition=reaction_condition,
                        dG_prime=dG_prime,
                        dG_std_dev=dG_std_dev
                    )
                    db.session.add(new_reaction)
                    db.session.commit()

                except Exception as e:
                    error_message = str(e)
                    logging.error(f"Error during prediction: {e}")

    return render_template(
        'prediction.html',
        form=form,
        dG_prime=dG_prime,
        dG_std_dev=dG_std_dev,
        error_message=error_message,
        conditions=conditions  # 传递条件字典以在模板中使用
    )

@app.route('/docs')
def docs():
    return render_template('docs.html')