import re

def extract_coefficients(expr, num_vars):
    coeffs = [0] * num_vars
    terms = re.findall(r'([+-]?\s*\d*)x(\d+)', expr)
    for coeff_str, var_idx_str in terms:
        var_idx = int(var_idx_str) - 1
        coeff_str = coeff_str.replace(' ', '')
        if coeff_str in ('', '+'):
            coeff = 1
        elif coeff_str == '-':
            coeff = -1
        else:
            coeff = int(coeff_str)
        coeffs[var_idx] = coeff
    return coeffs

def parse_lp_problem(file_path):
    """
    Usage:
    obj = parse_lp_problem(file_path)
    mode = obj['objective']['type']
    coefficients = obj['objective']['coefficients']

    constraints = obj['constraints']['lhs']
    signs = obj['constraints']['operators']
    rhs = obj['constraints']['rhs']
    """
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    variable_line = lines[-1]
    variable_indices = re.findall(r'x(\d+)', variable_line)
    num_vars = max(map(int, variable_indices))

    objective_line = lines[0].lower()
    is_min = 'min' in objective_line
    objective_expr = objective_line.split('=')[1]
    objective_coeffs = extract_coefficients(objective_expr, num_vars)
    
    lhs_coeffs = []
    rhs_values = []
    operators = []

    for line in lines[2:-1]:
        match = re.search(r'(<=|>=|=)', line)
        if not match:
            continue
        op = match.group()
        left, right = line.split(op)
        coeffs = extract_coefficients(left, num_vars)
        rhs = int(re.search(r'-?\d+', right.strip()).group())
        
        lhs_coeffs.append(coeffs)
        operators.append(op)
        rhs_values.append(rhs)
    
    obj =  {
        'objective': {
            'type': 'min' if is_min else 'max',
            'coefficients': objective_coeffs
        },
        'constraints': {
            'lhs': lhs_coeffs,
            'operators': operators,
            'rhs': rhs_values
        }
    }

    return fix_cst(obj)

def fix_cst(obj):
    for i in range(len(obj['constraints']['operators'])):

        if obj['constraints']['operators'][i] == '>=':

            obj['constraints']['lhs'][i] = [-coef for coef in obj['constraints']['lhs'][i]]
            obj['constraints']['rhs'][i] = -obj['constraints']['rhs'][i]
            obj['constraints']['operators'][i] = '<='
    
    return obj
