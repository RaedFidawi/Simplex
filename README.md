# Simplex solver for linear programs
Function is initially a minimization problem but it can change to max with the mode parameter. 

## Dependencies
- The python package Manim which is responsible for video generation should be installed.
Install using: 
`pip install manim`

- LaTeX Dependencies: 
texlive-latex 
texlive-collection-fontsrecommended 
dvipng 
dvisvgm
texlive-standalone
texlive-preview

Kindly install te requirements from requirements.txt by running:
pip install -r requirements.txt

## How to run

- change mode to min or max depending on your use case

- Add your constraints, which are of the form (expression <= value), RHS, and  objective function

- run the main.py using

`python main.py`