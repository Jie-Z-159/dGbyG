from flask import render_template, request, jsonify, send_file, current_app
from dgbygapp.forms import ReactionForm
from dgbygapp.utils import calculate_dg, process_batch_file, parse_compound, get_reaction_conditions
from dGbyG.api import Reaction, Compound
import os
import tempfile

def init_routes(app):
    @app.route('/')
    def home():
        return render_template('home.html')

    @app.route('/metabolite')
    def metabolite():
        return render_template('metabolite.html')

    @app.route('/prediction', methods=['GET', 'POST'])
    def prediction():
        form = ReactionForm()
        if form.validate_on_submit():
            # 处理单个反应的计算
            reaction = form.reaction.data
            dg = calculate_dg(reaction)
            return render_template('prediction.html', form=form, dg=dg, reaction=reaction)
        return render_template('prediction.html', form=form)

    @app.route('/parse_equation_route', methods=['POST'])
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

    @app.route('/calculate', methods=['POST'])
    def calculate():
        try:
            reaction_condition = request.form.get('reaction_condition')
            equation = request.form.get('equation')
            identifier = request.form.get('identifier_type')

            if not equation or not identifier or not reaction_condition:
                return jsonify({'error': 'Missing required parameters'}), 400

            conditions = get_reaction_conditions()

            if reaction_condition == 'custom':
                # 检查是否有全局自定义条件参数
                if all(key in request.form for key in ['global_pH', 'global_I', 'global_pMg', 'global_e_potential']):
                    # 使用全局自定义条件
                    custom_condition = {
                        'pH': float(request.form.get('global_pH')),
                        'I': float(request.form.get('global_I')),
                        'pMg': float(request.form.get('global_pMg')),
                        'e_potential': float(request.form.get('global_e_potential'))
                    }
                    
                    reaction = Reaction(equation, mol_type=identifier)
                    for compound, coeff in reaction.reaction.items():
                        compound.condition.update(custom_condition)
                    
                    dG_prime, dG_std_dev = reaction.transformed_standard_dGr_prime
                    dG_prime = float(dG_prime)
                    dG_std_dev = float(dG_std_dev)
                else:
                    # 使用原有的单个化合物自定义条件逻辑
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
                chosen_condition = conditions.get(reaction_condition)
                if not chosen_condition:
                    return jsonify({'error': 'Invalid reaction condition'}), 400

                reaction = Reaction(equation, mol_type=identifier)

                for compound, coeff in reaction.reaction.items():
                    compound.condition['pH'] = chosen_condition['pH']
                    compound.condition['I'] = chosen_condition['I']
                    compound.condition['pMg'] = chosen_condition['pMg']
                    compound.condition['e_potential'] = chosen_condition['e_potential']

                dG_prime, dG_std_dev = reaction.transformed_standard_dGr_prime
                dG_prime = float(dG_prime)
                dG_std_dev = float(dG_std_dev)
            
            return jsonify({
                'dG_prime': dG_prime, 
                'dG_std_dev': dG_std_dev,
                'message': 'Calculation completed successfully'
            })
            
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/batch_calculate', methods=['POST'])
    def batch_calculate():
        try:
            equations = request.form.get('equations', '').strip().split('\n')
            identifier = request.form.get('identifier_type')
            reaction_condition = request.form.get('reaction_condition')
            
            if not equations or not identifier or not reaction_condition:
                return jsonify({'error': 'Missing required parameters'}), 400

            # 获取全局自定义条件
            custom_condition = {
                'pH': float(request.form.get('global_pH', 7.0)),
                'I': float(request.form.get('global_I', 0.1)),
                'pMg': float(request.form.get('global_pMg', 7.0)),
                'e_potential': float(request.form.get('global_e_potential', 0.0))
            }

            results = []
            for equation in equations:
                equation = equation.strip()
                if not equation:
                    continue

                try:
                    # 创建反应对象并设置条件
                    reaction = Reaction(equation, mol_type=identifier)
                    for compound, coeff in reaction.reaction.items():
                        compound.condition.update(custom_condition)

                    # 计算ΔG和标准差
                    dG_prime, dG_std_dev = reaction.transformed_standard_dGr_prime
                    dG_prime = float(dG_prime)
                    dG_std_dev = float(dG_std_dev)

                    results.append({
                        'equation': equation,
                        'dG_prime': dG_prime,
                        'dG_std_dev': dG_std_dev,
                        'status': 'success'
                    })
                except Exception as e:
                    results.append({
                        'equation': equation,
                        'error': str(e),
                        'status': 'error'
                    })

            return jsonify({
                'results': results,
                'message': 'Batch calculation completed'
            })

        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api', methods=['GET', 'POST'])
    def api():
        if request.method == 'POST':
            if 'file' in request.files:
                file = request.files['file']
                if file.filename == '':
                    return jsonify({'error': 'No file selected'}), 400
                
                # 处理批量计算
                results = process_batch_file(file)
                return jsonify(results)
            else:
                # 处理单个反应的计算
                data = request.get_json()
                if not data or 'reaction' not in data:
                    return jsonify({'error': 'No reaction provided'}), 400
                
                reaction = data['reaction']
                dg = calculate_dg(reaction)
                return jsonify({'reaction': reaction, 'dg': dg})
        
        return render_template('api.html')

    @app.route('/citation')
    def citation():
        return render_template('citation.html')

    @app.route('/contact')
    def contact():
        return render_template('contact.html')

    @app.route('/faq')
    def faq():
        return render_template('faq.html') 