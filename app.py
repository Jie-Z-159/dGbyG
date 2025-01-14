import os
import sys
import click
import logging
from flask import Flask, render_template, request
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
app.secret_key = 'your_secret_key_for_session'  # WTF表单密钥
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
@app.cli.command()
def initdb():
    db.create_all()
    click.echo('Initialized database.')

# 首页和预测路由
@app.route('/', methods=['GET', 'POST'])
@app.route('/prediction', methods=['GET', 'POST'])
def prediction():
    form = ReactionForm()
    dG_prime = None
    dG_std_dev = None
    error_message = None

    if form.validate_on_submit():  # 检查表单是否提交且有效
        identifier = form.identifier_type.data
        equation = form.equation.data
        reaction_condition = form.reaction_condition.data
        
        try:
            # 创建 Reaction 对象，进行预测
            reaction = Reaction(equation, cid_type=identifier)
            dG_prime, dG_std_dev = reaction.transformed_standard_dGr_prime

            # 打印预测结果到终端（用于调试）
            print(f"Predicted ΔG′: {dG_prime}, Std Dev: {dG_std_dev}")

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
            # 捕获异常并设置错误消息
            error_message = str(e)
            dG_prime = None
            dG_std_dev = None
            logging.error(f"Error during prediction: {e}")

    return render_template(
        'prediction.html',
        form=form,
        dG_prime=dG_prime,
        dG_std_dev=dG_std_dev,
        error_message=error_message
    )

@app.route('/docs')
def docs():
    return render_template('docs.html')

if __name__ == '__main__':
    app.run(debug=True)
