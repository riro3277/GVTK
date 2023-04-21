import os
import math
import shutil
from datetime import date

class Testing():
    def __init__(self, CylDir, StandardDir):
        print("Im Right Heere", CylDir, StandardDir)
        self.StandardRoot = StandardDir
        self.RootDir = CylDir

    def Compare(self):

        StandardResults = self.StandardRoot + '/Lagrangian/StandardFiles/StandardFile0.txt'
        TestResults = self.RootDir + 'TestFiles/TestFile0.txt'

        StandardFile = open(StandardResults,'r')
        TestFile = open(TestResults,'r')

        StandardFile_lines = StandardFile.readlines()
        TestFile_lines = TestFile.readlines()

        self.DragErr = 0
        self.VelErr = 0
        self.ShearGradErr = 0
        self.BodyForceErr = 0
        self.PositionErr =  0
        self.TotalError = 0



        for i in range(len(TestFile_lines)):
            STDvals = StandardFile_lines[i].split(" ")
            STDtime = STDvals[0].split(":")[1]
            STDforce = STDvals[1].split(":")
            STDType = STDforce[0]
            STDvals = STDforce[1]
            STDVEC = STDvals.split(",")
            STDVEC = [float(STDVEC[0]), float(STDVEC[1]), float(STDVEC[2])]
            STDMag = math.sqrt((STDVEC[0]**2) + (STDVEC[1]**2) + (STDVEC[2]**2))

            Tvals = TestFile_lines[i].split(" ")
            Ttime = Tvals[0].split(":")[1]
            Tforce = Tvals[1].split(":")
            TType = Tforce[0]
            Tvals = Tforce[1]
            TVEC = Tvals.split(",")
            TVEC = [float(TVEC[0]), float(TVEC[1]), float(TVEC[2])]
            TMag = math.sqrt((TVEC[0]**2) + (TVEC[1]**2) + (TVEC[2]**2))

            error = abs(TMag - STDMag)

            self.TotalError += error
            if STDType == 'Drag' and STDType == TType:
                self.DragErr += error

            elif STDType == 'Velocity' and STDType == TType:
                self.VelErr += error

            elif STDType == 'ShearGradLift' and STDType == TType:
                self.ShearGradErr += error

            elif STDType == 'BodyForce' and STDType == TType:
                self.BodyForceErr += error

            elif STDType == 'Position' and STDType == TType:
                self.PositionErr += error


        if self.TotalError > 1e-3:
            self.Status = "FAILED"
        else:
            self.Status = "PASSED"


    def Documentation(self):
        TestDoc = self.RootDir + 'RegressionSummary.tex'
        TempDoc = self.StandardRoot + '/Lagrangian/TestingReportTemplate.tex'
        today = date.today()
        d1 = today.strftime("%m/%d/%Y")
        print("Regression Test Status", self.Status)
        shutil.copy(TempDoc, TestDoc)
        fin = open(TempDoc, "rt")

        fout = open(TestDoc, "w")

        for line in fin:


            if line.startswith('Date:'):
                fout.write(line.replace('0', d1))
            elif line.startswith('Velocity Error:'):
                fout.write(line.replace('0', str(self.VelErr)))
            elif line.startswith('Position Error:'):
                fout.write(line.replace('0', str(self.PositionErr)))
            elif line.startswith('Drag Error:'):
                fout.write(line.replace('0', str(self.DragErr)))
            elif line.startswith('ShearGradLift Error:'):
                fout.write(line.replace('0', str(self.ShearGradErr)))
            elif line.startswith('Body Force Error:'):
                fout.write(line.replace('0', str(self.BodyForceErr)))
            elif 'Regression Test Status:' in line:
                fout.write(line.replace('0', str(self.Status)))
            elif 'Summary/Actions:' in line:
                if self.Status == 'PASSED':
                    fout.write(line.replace('0', ' No action needed, regression test passed and it is safe to push new version to GitLab'))
                if self.Status == 'FAILED':
                    fout.write(line.replace('0', ' Regression test has failed meaning there are drastic changes in the results of the code. \nIf there is supposed to be drastic changes, please make sure you have conferred with the team that this is correct and then update the Standard File in gitlab by running the input-tracer-analytical-standard.dat for the cylindrical flow dataset.\nRun the regression test again with the new standard files to make sure it passes.'))

            else:
                fout.write(line)
        shutil.copy(TestDoc, self.StandardRoot + '/Docs/RegressionSummary.tex')




if __name__ == "__main__":
    baseDir = os.path.dirname(os.path.abspath(__file__))
    TestDir = baseDir + "/Cases/cylindrical_flow/"
    StandDir = baseDir + "/"
    print(baseDir)
    case = Testing(TestDir, StandDir)
    case.Compare()
    case.Documentation()
