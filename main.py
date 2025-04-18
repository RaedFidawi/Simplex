from manim import *
from fractions import Fraction
import copy

EPS = 1e-4

class Position:
    def __init__(self, column, row):
        self.column = column
        self.row = row

class NumberTypeclass:
    def zero(self): return 0
    def one(self): return 1
    def positive(self,x): return x > 0
    def iszero(self,x): return x == 0
    def nonnegative(self,x): return self.positive(x) or self.iszero(x)
    def coerce(self, x): return x
    def coerce_vec(self, x): return [self.coerce(xi) for xi in x]
    def coerce_mtx(self, x): return [self.coerce_vec(xi) for xi in x]
    

class RealFiniteTolerance(NumberTypeclass):
    def __init__(self, eps=1e-6):
        super(RealFiniteTolerance, self).__init__()
        self.eps = eps
        assert eps >= 0
        
    def zero(self): return 0.0
    def one(self): return 1.0
    def iszero(self,x): return abs(x) < self.eps
    def coerce(self, x): return float(x)

class RationalNumbers(NumberTypeclass):    
    def __init__(self):
        super(RationalNumbers, self).__init__()
        self._one = Fraction(1)
        self._zero = Fraction(0)
    def one(self): return self._one
    def zero(self): return self._zero
    def nonnegative(self,x): return x >= 0
    def coerce(self, x): return Fraction(x)

def _subtract_scaled_row(row1, row2, k, numclass):
    """row1 -= k*row2"""
    if numclass.iszero(k): return
    for i, row2_i in enumerate(row2):
        row1[i] -= k*row2_i

RESOLUTION_NO = "no"
RESOLUTION_SOLVED = "solved"
RESOLUTION_UNBOUNDED = "unbounded"
RESOLUTION_INCOMPATIBLE = "incompatible"

def epsilon_greater_than(a, b):
    return ((a > b) and not isclose(a, b))

def epsilon_greater_than_equal_to(a, b):
    return ((a > b) or isclose(a, b))

def epsilon_less_than(a, b):
    return ((a < b) and not isclose(a, b))

def epsilon_less_than_equal_to(a, b):
    return ((a < b) or isclose(a, b))

def isclose(a, b):
    return abs(a-b) <= EPS

class SimplexSolver:
    def __init__(self, a, b, c, basis=None, numclass=None, clean_c_row=False):
        assert len(a) == len(b) #number of rows of A must be the same as number of elements of b
        for aj in a:
            assert len(aj) == len(c) #each row of A must have the same len as c
            
        self.numclass = numclass or RealFiniteTolerance()
        self.a = a
        self.b = b
        self.c = c
        self.n = len(c)
        self.m = len(b)
        self.resolution = RESOLUTION_NO
        
        # If basis is not provided, initialize with None values
        if basis is None:
            self.basis = [None] * self.m
        else:
            self.basis = basis
            
        if clean_c_row: self._diagonalize_c_row()
        
        # For visualization, track history of states
        self.history = []
        self._save_state()
        
        # Additional variables for sm.py logic
        self.slack_rows = list(range(self.n, self.n + self.m))
        self.phase_one_optimization = False
        self.phase_one_row = None
        self.tableau = None
        self.last_pivot = None
        
    def _save_state(self):
        """Save the current state to history for visualization"""
        state = {
            'a': copy.deepcopy(self.a),
            'b': copy.deepcopy(self.b),
            'c': copy.deepcopy(self.c),
            'basis': copy.deepcopy(self.basis),
            'resolution': self.resolution
        }
        self.history.append(state)
        
    def _diagonalize_c_row(self):
        #diagonalize it (it has nonzero values now)
        c = self.c
        for j, i in enumerate(self.basis):
            if i is not None and not self.numclass.iszero(c[i]):
                _subtract_scaled_row(c, self.a[j], c[i], self.numclass)
                assert self.numclass.iszero(c[i])
                c[i] = self.numclass.zero()
        
    def vertex(self):
        v = [self.numclass.zero()] * self.n
        for j, i in enumerate(self.basis):
            if i is not None:  # Make sure the basis index is valid
                v[i] = self.b[j]
        return v
    
    def create_tableau(self, phase_one_optimization):
        tableau = []
        phase_one_row = [0] * (len(self.c) + self.m + 2)
        
        for i in range(self.m):
            if phase_one_optimization and self.b[i] < 0:
                slack_variables = [0] * self.m
                slack_variables[i] = -1.0
                tableau_row = [-1*x for x in self.a[i]] + slack_variables + [-1 * self.b[i]]
                tableau.append(tableau_row)
                phase_one_row = [a + b for a, b in zip(phase_one_row, tableau_row)]
            else:
                slack_variables = [0] * self.m
                slack_variables[i] = 1.0
                tableau_row = self.a[i] + slack_variables + [self.b[i]]
                tableau.append(tableau_row)
                
        final_row = [-1*x for x in self.c] + [0] * self.m + [0]
        tableau.append(final_row)
        
        return tableau, phase_one_row
    
    def select_pivot_element(self, tableau, phase_one_optimization, phase_one_row):
        pivot_element = Position(0, 0)
        no_solution = False
        
        if phase_one_optimization:
            pivot_element.column = phase_one_row.index(max(phase_one_row[:-1]))
        else:
            pivot_element.column = tableau[-1][:-1].index(min(tableau[-1][:-1]))
            
        ratios = []
        if pivot_element.column is not None:
            for r in range(len(tableau)-1):
                if tableau[r][pivot_element.column] > 0:
                    ratios.append(abs(tableau[r][-1] / tableau[r][pivot_element.column]))
                else:
                    ratios.append(float("inf"))
                    
            if all(i == float("inf") for i in ratios):
                no_solution = True
            else:
                row_min = min(ratios)
                row_min_indicies = [i for i,x in enumerate(ratios) if x == row_min]
                
                if (len(row_min_indicies) > 1):
                    least_variable = []
                    for j in row_min_indicies:
                        least_variable.append(self.slack_rows[j])
                    pivot_element.row = self.slack_rows.index(min(least_variable))
                else:
                    pivot_element.row = row_min_indicies[0]
        else:
            no_solution = True
            
        return no_solution, pivot_element
    
    def process_pivot_element(self, tableau, pivot_element, phase_one_optimization, phase_one_row):
        pri_mult = tableau[pivot_element.row][pivot_element.column]
        tableau[pivot_element.row] = [n / pri_mult for n in tableau[pivot_element.row]]
        tableau[pivot_element.row][pivot_element.column] = 1.0
        
        for i in range(len(tableau)):
            if i != pivot_element.row:
                sec_mult = tableau[i][pivot_element.column]
                pri_row = [j * sec_mult for j in tableau[pivot_element.row]]
                tableau[i] = [a - b for a, b in zip(tableau[i], pri_row)]
                tableau[i][pivot_element.column] = 0
                
        if phase_one_optimization:
            sec_mult = phase_one_row[pivot_element.column]
            pri_row = [j * sec_mult for j in tableau[pivot_element.row]]
            phase_one_row = [a - b for a, b in zip(phase_one_row, pri_row)]
            phase_one_row[pivot_element.column] = 0
            
        return tableau, phase_one_row
    
    def determine_answer(self, tableau):
        ans = [0] * self.n
        for i in range(self.n + self.m):
            if i < self.n and i in self.slack_rows:
                index = self.slack_rows.index(i)
                ans[i] = tableau[index][-1]
            elif i not in self.slack_rows and tableau[-1][i] == 0:
                for j in range(self.m-1):
                    if tableau[j][i] > 0:
                        return [-1]
            elif i < self.n:
                ans[i] = 0
        return ans
    
    def valid_answer(self, ans):
        invalid_answer = False
        for i in range(self.m):
            valid_ans = 0
            for j in range(self.n):
                valid_ans += self.a[i][j] * ans[j]
            if epsilon_greater_than(valid_ans, self.b[i]):
                invalid_answer = True
        if not all(epsilon_greater_than_equal_to(i, 0) for i in ans):
            invalid_answer = True
        return invalid_answer
    
    def step(self):
        # First time setup
        if self.tableau is None:
            # Check if we need phase 1 optimization (negative RHS)
            if any(bi < 0 for bi in self.b):
                self.phase_one_optimization = True
            self.tableau, self.phase_one_row = self.create_tableau(self.phase_one_optimization)
            
        # Check if we need to continue optimizing
        if (self.phase_one_optimization or 
            not all(epsilon_greater_than_equal_to(i, 0) for i in self.tableau[-1][:-1])):
            
            # Check if phase 1 is complete
            if self.phase_one_optimization and all(epsilon_less_than_equal_to(k, 0) for k in self.phase_one_row[:-1]):
                self.phase_one_optimization = False
                self.phase_one_complete = True
                # Check if we can proceed to phase 2
                if all(epsilon_greater_than_equal_to(i, 0) for i in self.tableau[-1][:-1]):
                    self.resolution = RESOLUTION_SOLVED
                    self._save_state()
                    return False
                    
            # Select pivot element
            no_solution, pivot_element = self.select_pivot_element(
                self.tableau, self.phase_one_optimization, self.phase_one_row
            )
            
            if no_solution:
                if hasattr(self, 'phase_one_complete') and self.phase_one_complete:
                    self.resolution = RESOLUTION_INCOMPATIBLE
                else:
                    self.resolution = RESOLUTION_UNBOUNDED
                self._save_state()
                return False
                
            # Update basis and perform pivot
            self.slack_rows[pivot_element.row] = pivot_element.column
            self.last_pivot = (pivot_element.row, pivot_element.column)
            
            self.tableau, self.phase_one_row = self.process_pivot_element(
                self.tableau, pivot_element, 
                self.phase_one_optimization, self.phase_one_row
            )
            
            # Update the state for visualization
            self._update_state_from_tableau()
            self._save_state()
            
            return True
        else:
            self.resolution = RESOLUTION_SOLVED
            self._save_state()
            return False

    def _update_state_from_tableau(self):
        """Update a, b, c from the current tableau state"""
        # Update a and b from the tableau (excluding slack variables)
        for j in range(self.m):
            self.a[j] = self.tableau[j][:self.n]
            self.b[j] = self.tableau[j][-1]
        
        # Update c from the last row of tableau (excluding slack variables)
        self.c = [-x for x in self.tableau[-1][:self.n]]
        
        # Update basis from slack_rows, properly mapping to original variables
        self.basis = [None] * self.m
        for row, col in enumerate(self.slack_rows):
            if row < self.m and col < self.n:  # Only map original variables
                self.basis[row] = col
            elif row < self.m and col >= self.n:  # Slack variable in basis
                # Find which constraint this slack variable corresponds to
                slack_idx = col - self.n
                if slack_idx < self.m:
                    self.basis[row] = None  # Or handle slack variables appropriately

    def solve(self):
        """Run the solver to completion"""
        while self.step():
            pass

# Manim visualization class
class SimplexVisualization(Scene):
    def construct(self):
        self.simplex_animation()
        
    def format_number(self, number):
        """Format numbers for display - handle fractions and floats"""
        if isinstance(number, Fraction):
            if number.denominator == 1:
                return str(number.numerator)
            return f"{number.numerator}/{number.denominator}"
        elif isinstance(number, float):
            if abs(number - round(number)) < 1e-10:
                return str(int(round(number)))
            return f"{number:.2f}"
        return str(number)
    
    def create_simplex_table(self, state, title=None, highlight_pivot=None):
        """Create a visualization of the simplex table for a given state"""
        a = state['a']
        b = state['b']
        c = state['c']
        basis = state['basis']
        
        m = len(b)
        n = len(c)
        
        table = VGroup()
        
        if title:
            title_text = Text(title, font_size=24)
            title_text.to_edge(UP)
            table.add(title_text)
        
        header = VGroup()
        corner_cell = Rectangle(height=0.6, width=0.8, color=WHITE)
        corner_text = Text("", font_size=16)
        corner_text.move_to(corner_cell.get_center())
        header.add(VGroup(corner_cell, corner_text))
        
        b_header = Rectangle(height=0.6, width=0.8, color=WHITE)
        b_text = Text("b", font_size=16)
        b_text.move_to(b_header.get_center())
        header.add(VGroup(b_header, b_text))
        
        for i in range(n):
            var_cell = Rectangle(height=0.6, width=0.8, color=WHITE)
            var_text = Text(f"x{i}", font_size=16)
            var_text.move_to(var_cell.get_center())
            header.add(VGroup(var_cell, var_text))
        
        header.arrange(RIGHT, buff=0.1)
        table.add(header)
        
        c_row = VGroup()
        c_label_cell = Rectangle(height=0.6, width=0.8, color=WHITE)
        c_label_text = Text("c", font_size=16)
        c_label_text.move_to(c_label_cell.get_center())
        c_row.add(VGroup(c_label_cell, c_label_text))
        
        empty_cell = Rectangle(height=0.6, width=0.8, color=WHITE)
        empty_text = Text("", font_size=16)
        empty_text.move_to(empty_cell.get_center())
        c_row.add(VGroup(empty_cell, empty_text))
        
        for i in range(n):
            cell = Rectangle(height=0.6, width=0.8, color=WHITE)
            text = Text(self.format_number(c[i]), font_size=16)
            text.move_to(cell.get_center())
            if c[i] < 0:
                cell.set_fill(RED_A, opacity=0.3)
            c_row.add(VGroup(cell, text))
        
        c_row.arrange(RIGHT, buff=0.1)
        table.add(c_row)
        
        for j in range(m):
            row = VGroup()
            basis_cell = Rectangle(height=0.6, width=0.8, color=WHITE)
            basis_text = Text(f"x{basis[j]}" if basis[j] is not None else "None", font_size=16)
            basis_text.move_to(basis_cell.get_center())
            row.add(VGroup(basis_cell, basis_text))
            
            b_cell = Rectangle(height=0.6, width=0.8, color=WHITE)
            b_text = Text(self.format_number(b[j]), font_size=16)
            b_text.move_to(b_cell.get_center())
            row.add(VGroup(b_cell, b_text))
            
            for i in range(n):
                cell = Rectangle(height=0.6, width=0.8, color=WHITE)
                text = Text(self.format_number(a[j][i]), font_size=16)
                text.move_to(cell.get_center())
                if highlight_pivot and highlight_pivot == (j, i):
                    cell.set_fill(YELLOW, opacity=0.5)
                if basis[j] == i:
                    cell.set_fill(GREEN_A, opacity=0.3)
                row.add(VGroup(cell, text))
            
            row.arrange(RIGHT, buff=0.1)
            table.add(row)
        
        table.arrange(DOWN, buff=0.1)
        return table
    
    def create_resolution_text(self, state):
        resolution = state['resolution']
        if resolution == RESOLUTION_SOLVED:
            return Text("Status: SOLVED", color=GREEN, font_size=24)
        elif resolution == RESOLUTION_UNBOUNDED:
            return Text("Status: UNBOUNDED", color=RED, font_size=24)
        elif resolution == RESOLUTION_INCOMPATIBLE:
            return Text("Status: INCOMPATIBLE", color=RED, font_size=24)
        else:
            return Text("Status: IN PROGRESS", color=YELLOW, font_size=24)

def visualize_simplex(a, b, c, basis=None, num=None, output_file="simplex_visualization"):
    num = num or RationalNumbers()
    
    class CustomSimplexVisualization(SimplexVisualization):
        def simplex_animation(self):
            solver = SimplexSolver(
                num.coerce_mtx(a), 
                num.coerce_vec(b), 
                num.coerce_vec(c),
                basis=basis,
                numclass=num
            )
            
            while solver.resolution == RESOLUTION_NO:
                solver.step()
            
            title = Text("Simplex Method Visualization", font_size=30)
            self.play(Write(title))
            self.wait(1)
            self.play(title.animate.to_edge(UP))
            
            # Helper function to format numbers for LaTeX
            def latex_format(number):
                if isinstance(number, Fraction):
                    if number.denominator == 1:
                        return f"{number.numerator}"
                    else:
                        return rf"\frac{{{number.numerator}}}{{{number.denominator}}}"
                elif isinstance(number, float):
                    if abs(number - round(number)) < 1e-9:
                        return f"{int(round(number))}"
                    else:
                        return f"{number:.2f}"
                else:
                    return f"{number}"
            
            # Build objective function string
            obj_terms = []
            for i in range(len(c)):
                coeff = c[i]
                term = f"{latex_format(coeff)}x_{{{i+1}}}"
                obj_terms.append(term)
            obj_str = " + ".join(obj_terms).replace("+ -", " - ")
            
            # Build constraint strings
            constraint_strs = []
            for i in range(len(a)):
                terms = []
                for j in range(len(a[i])):
                    coeff = a[i][j]
                    term = f"{latex_format(coeff)}x_{{{j+1}}}"
                    terms.append(term)
                lhs = " + ".join(terms).replace("+ -", " - ")
                rhs = latex_format(b[i])
                constraint_strs.append(f"{lhs} \\leq {rhs}")
            
            # Create LaTeX aligned environment
            lines = [
                r"\text{Maximize } & z = " + obj_str + r" \\",
                r"\text{Subject to:} & " + constraint_strs[0] + r" \\"
            ]
            for cs in constraint_strs[1:]:
                lines.append(r"& " + cs + r" \\")
            lines.append(r"& x_i \geq 0 \quad \forall i")
            
            full_tex = r"""
            \begin{aligned}
            %s
            \end{aligned}
            """ % " \\ ".join(lines)
            
            problem_tex = MathTex(full_tex, font_size=24)
            problem_tex.next_to(title, DOWN, buff=0.5)
            
            self.play(Write(problem_tex))
            self.wait(2)
            self.play(FadeOut(problem_tex))
            
            # Rest of the animation code remains the same...
            initial_table = self.create_simplex_table(solver.history[0], "Initial Simplex Table")
            self.play(FadeIn(initial_table))
            self.wait(2)
            
            current_table = initial_table
            for i in range(1, len(solver.history)):
                state = solver.history[i]
                
                pivot = None
                if hasattr(solver, 'last_pivot') and i > 0:
                    pivot = solver.last_pivot
                
                step_title = f"Step {i}"
                if state['resolution'] != RESOLUTION_NO:
                    step_title = f"Final State: {state['resolution'].capitalize()}"
                    
                new_table = self.create_simplex_table(state, step_title, highlight_pivot=pivot)
                
                explanation = None
                if pivot:
                    row, col = pivot
                    explanation = Text(
                        f"Pivot: row {row+1}, column {col+1} (x{col})",
                        font_size=18
                    )
                    explanation.next_to(new_table, DOWN, buff=0.5)
                
                self.play(Transform(current_table, new_table))
                
                if explanation:
                    self.play(Write(explanation))
                    self.wait(2)
                    self.play(FadeOut(explanation))
                else:
                    self.wait(2)
            
            resolution_text = self.create_resolution_text(solver.history[-1])
            resolution_text.next_to(current_table, DOWN, buff=0.5)
            self.play(Write(resolution_text))
            
            self.wait(2)
            
            if solver.history[-1]['resolution'] == RESOLUTION_SOLVED:
                solution = solver.vertex()
                solution_text = Text(
                    f"Optimal Solution: x = [{', '.join(self.format_number(x) for x in solution)}]",
                    font_size=20
                )
                solution_text.next_to(resolution_text, DOWN, buff=0.5)
                
                optimal_value = sum(c[i] * solution[i] for i in range(len(c)))
                optimal_text = Text(
                    f"Optimal Value: z = {self.format_number(optimal_value)}",
                    font_size=20
                )
                optimal_text.next_to(solution_text, DOWN, buff=0.3)
                
                self.play(Write(solution_text))
                self.play(Write(optimal_text))
            
            self.wait(3)
    
    config.output_file = f"{output_file}.mp4"
    scene = CustomSimplexVisualization()
    scene.render()

if __name__ == "__main__":

    mode ='min'

    a = [
        [-8, -16],  # First constraint (≤)
        [-60, -40],  # Second constraint (≥)
        [-2, -2]
    ]

    b = [-200, -960, -40]

    c = [60, 50]  # Minimize this objective function

    if mode == 'max':
        visualize_simplex(a, b, c, output_file="simplex_example1")
    else:
        c = [-1*i for i in c]
        visualize_simplex(a, b, c, output_file="simplex_example1")