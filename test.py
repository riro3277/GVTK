import os
#runs regression test without making test fles. this is only for after locally creatinng test files by using ./Lagrangian/executeMainGVTK.py in test mode
dir = os.path.dirname(os.path.abspath(__file__))
with open(dir + "/Lagrangian/testModule.py" ) as f:
    exec(f.read())
