#-------------------------------------------------------------------------------
# import modules
#-------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import sys, os, vtk
import matplotlib.pyplot as plt
import matplotlib.image as image
import moduleExposureRisk as M_EXP


import json

#-------------------------------------------------------------------------------
# calculate viral decay
#-------------------------------------------------------------------------------
def viabilityExpDecay(half_life, t):
    return np.exp(-t*np.log(2)/half_life)

#-------------------------------------------------------------------------------
# read .dat file
#-------------------------------------------------------------------------------
inputFile     = sys.argv[1]
exposureInput = M_EXP.ExposureInputs(inputFile)
exposureInput.readInputFile()

#-------------------------------------------------------------------------------
# construct path to files
#-------------------------------------------------------------------------------
results_dir   = exposureInput.m_ResultsDirectory
series_name   = '_'.join(os.listdir(results_dir)[0].split('.')[0].split('_')[:-1])+'_'
particle_ext  = '.vtk'
particle_name = results_dir + series_name

#-------------------------------------------------------------------------------
# subject data import
#-------------------------------------------------------------------------------
emit_subs = exposureInput.m_EmittingSubjects.split(';')
emit_subs = [x.strip(' ') for x in emit_subs]
emit_subs = np.fromiter(emit_subs, dtype=int)

exp_subs     = exposureInput.m_ExposedSubjects.split(';')
exp_subs     = [x.strip(' ') for x in exp_subs]
exp_subs     = np.fromiter(exp_subs, dtype=int)
num_exp_subs = len(exp_subs)

head_centers     = M_EXP.arrayParser(exposureInput.m_HeadCenters, 2, float)
exp_head_centers = head_centers[exp_subs]

#-------------------------------------------------------------------------------
# exposure zone setup
#-------------------------------------------------------------------------------
exposureZones = M_EXP.createZoneDictionary(exposureInput.m_ExposureZoneFolder, exposureInput.m_SubjectNames, exp_subs)
exposedNames  = list(exposureZones.keys())

# print(json.dumps(exposureZones, sort_keys=False, indent=4))
# print(exposureZones)
# sys.exit()

# num_exp_zones  = exposureInput.m_NumExposureZones
# zone_radii     = np.linspace(exposureInput.m_Zone0Radius, exposureInput.m_Zone0Radius+0.1, num_exp_zones+1)[1:]
# exp_zone_shape = exposureInput.m_ExposureZoneShape
# if exp_zone_shape == 'cone':
#     orientations = exposureInput.m_Orientations
#     zone_angle   = exposureInput.m_ZoneAngle
# elif exp_zone_shape == 'cylinder':
#     orientations = exposureInput.m_Orientations
#     zone_angle   = exposureInput.m_ZoneAngle
#     zone_height  = exposureInput.m_ZoneHeight

#-------------------------------------------------------------------------------
# 2D exposure grid setup
#-------------------------------------------------------------------------------
isGenerateHeatmaps = exposureInput.m_GenerateHeatmaps
if isGenerateHeatmaps:
    resolution     = exposureInput.m_GridResolution
    x_range        = exposureInput.m_GridXRange
    y_range        = exposureInput.m_GridYRange
    numXGridPoints = round((x_range[1]-x_range[0])/resolution)
    numYGridPoints = round((y_range[1]-y_range[0])/resolution)
    x_arr          = np.linspace(x_range[0], x_range[1], numXGridPoints+1)
    y_arr          = np.linspace(y_range[0], y_range[1], numYGridPoints+1)
    includeBackground  = exposureInput.m_IncludeBackground
    if includeBackground:
        backImg = image.imread(exposureInput.m_BackgroundImage)
    image_iter  = 0

#-------------------------------------------------------------------------------
# viability tracking set up
#-------------------------------------------------------------------------------
decay_rate   = exposureInput.m_DecayRate
viability_dt = exposureInput.m_ViabilityDelta
linear_var   = exposureInput.m_LinearViability

#-------------------------------------------------------------------------------
# get number of particles
#-------------------------------------------------------------------------------
reader      = vtk.vtkPolyDataReader()
reader.SetFileName(particle_name+'0'+particle_ext)
reader.Update()
num_points  = reader.GetOutput().GetNumberOfPoints()

#-------------------------------------------------------------------------------
# calculate exposure risk for each subject
#-------------------------------------------------------------------------------
exp_indices = np.zeros(num_exp_subs)
num_files   = len(os.listdir(results_dir))
file_step   = exposureInput.m_FileStepSize
time_steps  = np.linspace(0, (num_files-1)*file_step, num_files, dtype='int')
# time_steps  = np.linspace(0, 600, 601, dtype='int')

for time in time_steps:

    #-------------------------
    # print status indication
    #-------------------------
    # if int(time/file_step)+1 == num_files:
    #     print('Step', int(time/file_step)+1, 'of', num_files)
    # else:
    #     print('Step', int(time/file_step)+1, 'of', num_files, end='\r')

    #----------------------------------------
    # update particle file for loop progress
    #----------------------------------------
    reader.SetFileName(particle_name+str(time)+particle_ext)
    reader.Update()

    #--------------------------------------
    # generate array of zeros for heatmaps
    #--------------------------------------
    if isGenerateHeatmaps:
        exposureGrid = np.zeros((numXGridPoints, numXGridPoints))

    #------------------
    # update viability
    #------------------
    if linear_var:
        viability -= decay_rate*viability_dt
    else:
        viability  = viabilityExpDecay(decay_rate, time*viability_dt)

    for i in range(num_points):
        point = reader.GetOutput().GetPoint(i)

        if int(time/file_step)+1 == num_files:
            print('\x1b[2K\rStep', int(time/file_step)+1, 'of', num_files, '| Particle:', i)
        else:
            print('\x1b[2K\rStep', int(time/file_step)+1, 'of', num_files, '| Particle:', i, end='\r')

        #------------------------------------------
        # calculate exposure risk for each subject
        #------------------------------------------
        for j in range(num_exp_subs):
            # if exp_zone_shape == 'cone':
            #     exp_indices[j] += M_EXP.exposureRiskCone(point, exp_head_centers[j], zone_radii, zone_angle, num_exp_zones, viability, orientations[j])
            # elif exp_zone_shape == 'cylinder':
            #     exp_indices[j] += M_EXP.exposureRiskCyl(point, exp_head_centers[j], zone_radii, zone_angle, zone_height, num_exp_zones, viability, orientations[j])
            # elif exp_zone_shape == 'sphere':
            #     exp_indices[j] += M_EXP.exposureRiskSphere(point, exp_head_centers[j], zone_radii, num_exp_zones, viability)
            exp_indices[j] += M_EXP.calculateExposureRisk(exposureZones[exposedNames[j]], point, viability)

        #----------------------------------------------
        # find bin in exposureGrid that holds particle
        # increment bin count by 1
        #----------------------------------------------
        if isGenerateHeatmaps:
            x = np.searchsorted(x_arr, point[0], side='left')-1
            y = np.searchsorted(y_arr, point[1], side='left')-1
            if x == -1:
                x = 0
            if x == numXGridPoints:
                x -= 1
            if y == -1:
                y = 0
            if y == numYGridPoints:
                y -= 1
            exposureGrid[x,y] += 1

    #--------------------------------
    # plot exposureGrid as a heatmap
    #--------------------------------
    if isGenerateHeatmaps:
        if includeBackground:
            imFileName = exposureInput.m_AnimationFilename + str(image_iter) + '.png'
            fig, ax = plt.subplots()
            ax.imshow(backImg, extent=(0,numXGridPoints,0,numYGridPoints), aspect='auto')
            plt.pcolormesh(exposureGrid, cmap='viridis', edgecolors='w', linewidth=0.01, vmin=0, vmax=500, alpha=0.8)
        else:
            plt.figure()
            plt.pcolormesh(exposureGrid, cmap='viridis', edgecolors='w', linewidth=0.01, vmin=0, vmax=500)
        plt.title('Exposure Grid')
        plt.xticks(ticks=range(numXGridPoints+1), labels=np.round((x_arr), 2), rotation=90)
        plt.yticks(ticks=range(numYGridPoints+1), labels=np.round((y_arr), 2))
        clb = plt.colorbar(label='Number of Particles in Bin')
        clb.set_ticks([0, 100, 200, 300, 400, 500])
        clb.set_ticklabels(['0','100','200','300','400', '≥500'])
        plt.xlabel('m'), plt.ylabel('m')
        plt.savefig(imFileName)
        plt.close()
        image_iter += 1

#-------------------------------------------------------------------------------
# normalize exposure indices
#-------------------------------------------------------------------------------
exp_indices /= num_files * num_points

#-------------------------------------------------------------------------------
# output exposure risk calculations to stored file
#-------------------------------------------------------------------------------
mainFileName  = exposureInput.m_ExpFileName
headerNames   = exposureInput.m_HeaderNames
headerEntries = exposureInput.m_HeaderEntries
data = pd.DataFrame()

for i in range(len(headerNames)):
    data.loc[0, headerNames[i]] = headerEntries[i]

for i in range(len(emit_subs)+num_exp_subs):
    data.loc[0, 'Occupant '+str(i)] = np.nan

for i in range(num_exp_subs):
    data.loc[0, 'Occupant '+str(exp_subs[i])] = exp_indices[i]

if os.path.isfile(mainFileName):
    main = pd.read_csv(mainFileName, header=0)
    main = pd.concat([main, data], sort=False)
    main.to_csv(mainFileName, index=False)
else:
    data.to_csv(mainFileName, index=False)
