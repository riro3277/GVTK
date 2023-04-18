rm *.o
rm *.mod
rm *.so
gfortran -c -fPIC -O3 boundaryConditions.f90
gfortran -c -fPIC -O3 mymath.f90
gfortran -c -fPIC -O3 analyticalFlowLibraryFT.f90
gfortran -c -fPIC -O3 velocity.f90
gfortran -c -fPIC -O3 outputFormat.f90
ar crs libMesh.a mymath.o
gfortran -c -fPIC -O3 mesh.f90
ar crs libIntegration.a velocity.o boundaryConditions.o analyticalFlowLibraryFT.o mymath.o mesh.o
gfortran -c -fPIC -O3 integration.f90
gfortran -c -fPIC -O3 cgFields.f90
ar crs libMainCall.a integration.o cgFields.o mymath.o analyticalFlowLibraryFT.o boundaryConditions.o velocity.o mesh.o outputFormat.o
f2py3 --quiet -m mesh -h meshSgnFile.pyf mesh.f90 --overwrite-signature
f2py3 --quiet -m integration -h integrationSgnFile.pyf integration.f90 --overwrite-signature
f2py3 --quiet -m main -h mainSgnFile.pyf main.f90 --overwrite-signature
f2py3 --quiet --opt='-O3' -c --fcompiler=gnu95 --f90flags=-w meshSgnFile.pyf mesh.f90 -L. -lMesh
f2py3 --quiet --opt='-O3' -c --fcompiler=gnu95 --f90flags=-w integrationSgnFile.pyf integration.f90 -L. -lIntegration
f2py3 --quiet --opt='-O3' -c --fcompiler=gnu95 --f90flags=-w mainSgnFile.pyf main.f90 -L. -lMainCall -llapack -lrefblas -ltmglib
