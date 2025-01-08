from flask import Flask, render_template, request, jsonify
from dGbyG.api import Compound, Reaction
from dGbyG.utils.ChemFunc import parse_equation
from forms import ReactionForm

app = Flask(__name__)
app.secret_key = 'your_secret_key_for_session'#wtf表单密钥
# 首页和大子页路由
@app.route('/')
@app.route('/prediction',methods=['GET', 'POST'])
def prediction():
    form = ReactionForm()
    if form.validate_on_submit():  # 检查表单是否提交且有效
        identifier = form.identifier_type.data
        equation = form.equation.data
        cell_compartment = form.cell_compartment.data

        # 处理表单数据
        return f"Identifier Type: {identifier}, Equation: {equation}, Cell Compartment: {cell_compartment}"

    return render_template('prediction.html',form=form)

@app.route('/docs')
def docs():
    return render_template('docs.html')



