import os
import sys
import click
import logging
import re
from flask import Flask, render_template, request,jsonify
from flask_sqlalchemy import SQLAlchemy
from dGbyG.api import Compound, Reaction
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
    dG_prime = db.Column(db.Float, nullable=True)
    dG_std_dev = db.Column(db.Float, nullable=True) 

    def __repr__(self):
        return '<ReactionModel %r>' % self.equation

# 定义固定条件
default_T = 298.15
default_I = 0.25
default_pMg = 14.0


@app.cli.command('initdb')
def initdb():
    db.create_all()
    click.echo('Initialized database.')

def parse_compound(compound):
    # 用正则提取数字（如果有的话）和化合物名字
    match = re.match(r'^\s*(\d*)(.+?)\s*$', compound)
    if match:
        num = match.group(1)
        name = match.group(2).strip()
        num = int(num) if num else 1  # 没有数字就默认是1
        return num, name
    else:
        raise ValueError(f"Cannot parse compound: {compound}")

@app.route('/parse_equation', methods=['POST'])
def parse_equation_route():
    equation = request.form.get('equation', '').strip()
    if ' = ' not in equation:
        return jsonify({'error': 'Equation must contain " = " separator'}), 400
    
    lhs, rhs = equation.split(' = ', 1)
    substrates = [c.strip() for c in lhs.split('+')]
    products = [c.strip() for c in rhs.split('+')]

    compound_names = []
    stoichiometries = []

    for c in substrates:
        coeff, name = parse_compound(c)
        compound_names.append(name)
        stoichiometries.append(-coeff)  # 左边是消耗，负数

    for c in products:
        coeff, name = parse_compound(c)
        compound_names.append(name)
        stoichiometries.append(+coeff)  # 右边是生成，正数

    return jsonify({'compounds': compound_names, 'stoichiometries': stoichiometries})

@app.route('/', methods=['GET'])
def home():
    return render_template('home.html')

@app.route('/prediction', methods=['GET'])
def prediction():
    form = ReactionForm()
    return render_template('prediction.html',form = form)

@app.route('/calculate',methods = ['POST'])
def calculate():
    form = ReactionForm()
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

    if not form.validate_on_submit():
        return jsonify({'error': 'Form validation failed', 'errors': form.errors}), 400
        
    try:
        reaction_condition = form.reaction_condition.data
        equation = form.equation.data
        identifier = form.identifier_type.data
        
        if reaction_condition == 'custom':
            lhs, rhs = equation.split(' = ',1)
            substrates = [c.strip() for c in lhs.split('+')]
            products = [c.strip() for c in rhs.split('+')]

            compound_names = []
            stoichiometries = []

            for c in substrates:
                coeff, name = parse_compound(c)
                compound_names.append(name)
                stoichiometries.append(-coeff)

            for c in products:
                coeff, name = parse_compound(c)
                compound_names.append(name)
                stoichiometries.append(coeff)

            # 用索引区分每个化合物的条件
            custom_conditions = {}
            for index, compound in enumerate(compound_names):
                pH = float(request.form.get(f'custom_conditions-{index}-pH', 7.0))
                I = float(request.form.get(f'custom_conditions-{index}-I', 0.1))
                pMg = float(request.form.get(f'custom_conditions-{index}-pMg', 7.0))
                e_potential = float(request.form.get(f'custom_conditions-{index}-e_potential', 0.0))
                custom_conditions[index] = {
                    'pH': pH,
                    'I': I,
                    'pMg': pMg,
                    'e_potential': e_potential
                }

            # 创建化合物对象并设置自定义条件（用索引区分）
            compound_objs = []
            for index, compound_name in enumerate(compound_names):
                compound_obj = Compound(compound_name, input_type=identifier)
                compound_obj.condition = custom_conditions.get(index, {})
                compound_objs.append(compound_obj)

            equation_dict = dict(zip(compound_objs, stoichiometries))
            reaction = Reaction(equation_dict, mol_type='compound')
            dG_prime, dG_std_dev = reaction.transformed_standard_dGr_prime
            dG_prime = float(dG_prime)
            dG_std_dev = float(dG_std_dev)
        else:
            reaction_condition = form.reaction_condition.data
            identifier = form.identifier_type.data
            equation = form.equation.data

            chosen_condition = conditions.get(reaction_condition)
            reaction = Reaction(equation, mol_type=identifier)

            for compound,coeff in reaction.reaction.items():
                compound.condition['pH'] = chosen_condition['pH']
                compound.condition['I'] = chosen_condition['I']
                compound.condition['pMg'] = chosen_condition['pMg']
                compound.condition['e_potential'] = chosen_condition['e_potential']
            dG_prime, dG_std_dev = reaction.transformed_standard_dGr_prime
            dG_prime = float(dG_prime)
            dG_std_dev = float(dG_std_dev)

        # 保存结果到数据库
        reaction_model = ReactionModel(
            identifier=identifier,
            equation=equation,
            dG_prime=dG_prime,
            dG_std_dev=dG_std_dev
        )
        db.session.add(reaction_model)
        db.session.commit()
        
        return jsonify({
            'dG_prime': dG_prime, 
            'dG_std_dev': dG_std_dev,
            'message': 'Calculation completed successfully'
        })
        
    except Exception as e:
        db.session.rollback()
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, use_reloader=True)
