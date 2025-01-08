from flask import Flask, render_template, request, jsonify
from dGbyG.api import Compound, Reaction
from dGbyG.utils.ChemFunc import parse_equation

app = Flask(__name__)
app.secret_key = 'your_secret_key_for_session'#wtf表单密钥
# 首页和大子页路由
@app.route('/')
@app.route('/prediction')
def prediction():
    return render_template('prediction.html')

@app.route('/docs')
def docs():
    return render_template('docs.html')



