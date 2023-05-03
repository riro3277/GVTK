import os
import math
import shutil
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt


class Testing():
    def __init__(self, CylDir, StandardDir):
        #Standard dir is the base gvtktest folder location
        #Root dir is the location of the ./cylinder_flow folder
        self.StandardRoot = StandardDir
        self.RootDir = CylDir

    def Compare(self):
        #Reading standard files and test files from both locations for comparison
        StandardResults = self.StandardRoot + '/Lagrangian/StandardFiles/StandardFile0.txt'
        TestResults = self.RootDir + 'TestFiles/TestFile1.txt'

        StandardFile = open(StandardResults,'r')
        TestFile = open(TestResults,'r')

        StandardFile_lines = StandardFile.readlines()
        TestFile_lines = TestFile.readlines()
        #Initizalzing particle parameter error measurements
        self.DragErr = 0
        self.VelErr = 0
        self.ShearGradErr = 0
        self.BodyForceErr = 0
        self.PositionErr =  0
        self.TotalError = 0
        #Initial arrays for particle position chart that shows standard and test results
        xtest = []
        ytest =[]
        ztest = []
        xstd = []
        ystd = []
        zstd = []
        #Comparison loop that extracts each value from the test files and standard files for the error difference
        #STDvals are the values from the standard Files were each dimension is stored in STDVEC.
        #Tvals are the values are the values form the test file
        #Type values are the name of the parameter (i.e. Drag or position)
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
            #Computes error for each value
            error = abs(TMag - STDMag)
            #appends error to each specific paramter type to create a summed error across all the time steps
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
                #get posiiton vectors in each dimension for test and standard data for plotting
                xstd.append(STDVEC[0])
                ystd.append(STDVEC[1])
                zstd.append(STDVEC[2])
                xtest.append(TVEC[0])
                ytest.append(TVEC[1])
                ztest.append(TVEC[2])

        #if the total summed error across all parameters is non-zero, regression test is failed
        if self.TotalError > 1e-3:
            self.Status = "FAILED"
        else:
            self.Status = "PASSED"

        #Plot postition of test and standard results
        ax = plt.axes(projection='3d')
        ax.scatter3D(xstd, ystd, zstd, c = 'red', label = 'Standard Results')
        ax.scatter3D(xtest, ytest, ztest, c = 'green', label = 'Test Results')
        plt.legend()
        plt.title("Particle Position of Standard and Test Results")
        plt.savefig(self.StandardRoot + '/Docs/Position.png', dpi =300)


    def Documentation(self):
        #.tex file that the regression test being completed will inform
        TestDoc = self.RootDir + 'RegressionSummary.tex'
        #Template is edited to include information of the specific regression test
        TempDoc = self.StandardRoot + '/Lagrangian/TestingReportTemplate.tex'
        today = datetime.now()
        d1 = today.strftime("%m/%d/%Y %H:%M:%S")
        print("Regression Test Status:", self.Status)
        print("Please see summary within Docs folder for more information")
        #A failed regression test will return a divde by zero error for CI in GitLab in order to notify that the regression test has failed. Documentation of failure is still created
        if self.Status == "FAILED":
            print("Divide by zero error is for GitLab CI, supposed to be there")
        shutil.copyfile(TempDoc, TestDoc)
        fin = open(TempDoc, "rt")

        fout = open(TestDoc, "w")
        #editing template to include information about this specific regression test and then copying it to the ./Docs folder
        #Includes a summary of the test, the summed errors for each parameter, next steps if the test passes or fails, and the position plot of the standard data against the test data
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
        fin.close()
        fout.close()
        shutil.copyfile(TestDoc, self.StandardRoot + '/Docs/RegressionSummary.tex')
        #Divide by zero error for GitLab CI
        if self.Status == "FAILED":
            return(1/0)




if __name__ == "__main__":
    baseDir = os.path.dirname(os.path.abspath(__file__))
    TestDir = baseDir + "/Cases/cylindrical_flow/"
    StandDir = baseDir + "/"
    case = Testing(TestDir, StandDir)
    case.Compare()
    case.Documentation()
