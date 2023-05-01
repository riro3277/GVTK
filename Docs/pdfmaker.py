import os
dir = os.path.dirname(os.path.abspath(__file__))
file = dir + "/RegressionSummary.tex"
os.system('pdflatex ' + file)
