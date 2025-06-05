import re
from dGbyG.api import Compound, Reaction

# 定义固定条件
default_T = 298.15
default_I = 0.25
default_pMg = 14.0

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

def get_reaction_conditions():
    return {
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

def calculate_dg(reaction_str):
    """计算单个反应的ΔG值"""
    try:
        reaction = Reaction(reaction_str, mol_type='compound')
        dG_prime, dG_std_dev = reaction.transformed_standard_dGr_prime
        return {
            'dG_prime': float(dG_prime),
            'dG_std_dev': float(dG_std_dev)
        }
    except Exception as e:
        return {'error': str(e)}

def process_batch_file(file):
    """处理批量计算文件"""
    results = []
    for line in file:
        line = line.decode('utf-8').strip()
        if not line:
            continue
        try:
            result = calculate_dg(line)
            result['reaction'] = line
            results.append(result)
        except Exception as e:
            results.append({
                'reaction': line,
                'error': str(e)
            })
    return results 