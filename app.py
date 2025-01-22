import os
import sys
import click
import logging
import json
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

@app.route('/', methods=['GET'])
@app.route('/prediction', methods=['GET'])
def prediction():
    form = ReactionForm()
    return render_template('prediction.html',form = form)


@app.route('/calculate', methods=['POST'])
def calculate():
    form = ReactionForm()
    if form.validate_on_submit():
        reaction_condition = form.reaction_condition.data  # 获取反应条件
        equation = form.equation.data  # 获取化学方程式
        identifier = form.identifier_type.data  # 获取标识符

        # 处理非自定义条件
        chosen_condition = conditions.get(reaction_condition)
        reaction = Reaction(equation, cid_type=identifier)

        # 根据条件设置每个化合物的条件
        for compound, coeff in reaction.reaction.items():
            compound.condition['pH'] = chosen_condition['pH']
            compound.condition['I'] = chosen_condition['I']
            compound.condition['pMg'] = chosen_condition['pMg']
            compound.condition['e_potential'] = chosen_condition['e_potential']

        dG_prime, dG_std_dev = reaction.transformed_standard_dGr_prime
        dG_prime = float(dG_prime)
        dG_std_dev = float(dG_std_dev)

    return jsonify({'dG_prime': dG_prime, 'dG_std_dev': dG_std_dev})


@app.route('/docs')
def docs():
    return render_template('docs.html')