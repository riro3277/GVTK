import vtk, sys

def vtkTetrahedralizeFOAMData(a_FoamFile, a_TetFile, a_ReturnData=False, a_FoamBlockID=0):

    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(a_FoamFile)
    reader.SkipZeroTimeOn()
    reader.Update()
    foamData = reader.GetOutput()
    numBlocks = foamData.GetNumberOfBlocks()
    if numBlocks == 1:
        blockData = foamData.GetBlock(0)
    else:
        blockData = foamData.GetBlock(a_FoamBlockID)

    tetFilter = vtk.vtkDataSetTriangleFilter()
    tetFilter.SetInputData(blockData)
    tetFilter.TetrahedraOnlyOn()
    tetFilter.Update()

    tetData = tetFilter.GetOutput()

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(a_TetFile)
    writer.SetDataModeToBinary()
    writer.SetInputData(tetData)
    writer.Update()
    writer.Write()

    if a_ReturnData == True:
        return tetData

def vtkVerifyAllTets(a_TetData):

    numCells    = a_TetData.GetNumberOfCells()
    tet         = 0

    for cell in range(numCells):
        if a_TetData.GetCellType(cell) == vtk.VTK_TETRA:
            tet = tet + 1

    if tet == numCells:
        print('All cells are tetrahedral simplices')
    else:
        print('Grid has mixed cell type simplices')

if __name__=='__main__':

    inFile      = sys.argv[1].split('/') #path/to/flow/data/flow-data-file.foam
    filePath    = '/'.join(inFile[:-1])+'/'
    foamFile    = inFile[-1]
    tetFile     = foamFile.split('.')[0]+'.vtu'
    data        = vtkTetrahedralizeFOAMData(filePath+foamFile, filePath+tetFile, a_ReturnData=True)
    vtkVerifyAllTets(data)
