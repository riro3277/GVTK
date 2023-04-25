import os
dir = os.path.dirname(os.path.abspath(__file__))
gvtk = dir +"/Lagrangian/executeMainGVTK.py"
v = os.system("python3" + " " +  gvtk + " test")
