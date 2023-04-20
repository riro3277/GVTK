import os
print(os.path.dirname(os.path.abspath(__file__)))
dir = os.path.dirname(os.path.abspath(__file__))
with open(dir + "/Lagrangian/testModule.py" ) as f:
    exec(f.read())
