from paraview.simple import *

"""python script to postprocess Green's functions in paraview. To run this, open paraview, go to tools --> python shell --> run script and select this file.
Note that your data files from the SPECFEMX runs should be named with a convention that includes the subfault number in the .case file names so you can loop over them easily as shown below. 
You will need a file containing the locations of the points at which you want to save data -- generally the locations of GPS stations and InSAR data -- in .vtk format. Make sure that the locations in this file use the same coordinate system as the mesh used to calculate the GFs."""

numfaults = 180 # how many patches to process. change this to the number of subfaults you've used for your GFs

numproc = 40 #each run is assumed to have 40 processors. change this number otherwise
ranger = [str(x).zfill(2) for x in range(numproc)] #list of strings with processor numbers

#list of modeling results from all subfaults. we need both strikeslip and dipslip GFs
models1=["%d_ds" % x for x in range(0,numfaults)]
models2=["%d_ss" % x for x in range(0,numfaults)]
models=models1+models2

print models
for model in models:
    print model

    #read in all case files for this particular modeling result: subfault + strikeslip or dipslip
    case_files=[EnSightReader(CaseFileName='your_directory_path/your_file_name%s_proc%s.case' %(model,x)) for x in ranger]
    
    #we only want to read in the displacement result
    for cfile in case_files:
        cfile.PointArrays = ['displacement']

    #group and merge all processors
    groupDatasets1 = GroupDatasets(Input=case_files)
    mergeBlocks1 = MergeBlocks(Input=groupDatasets1)

    # extract the model surfaces
    extractSurface1 = ExtractSurface(Input=mergeBlocks1)

    # clip the top surface
    clip1 = Clip(Input=extractSurface1)
    clip1.ClipType = 'Plane'
    clip1.Scalars = [None, '']
    clip1.ClipType.Origin = [0.0, 0.0, -1.0] #this should be X center coordinate, Y center coordinate, and the Z coordinate that is just below the lowest point on the surface of the mesh. You may need to experiment with clipping manually to find this number.
    clip1.ClipType.Normal = [0.0, 0.0, 1.0]

    # flatten the topographic surface using the calculator
    calculator1 = Calculator(Input=clip1)
    calculator1.Function = ''
    calculator1.CoordinateResults = 1
    calculator1.Function = 'coords-(0*iHat+0*jHat+coordsZ*kHat)'

    # read in the station coordinates. we want Green's functions only at these points
    nepal_station_coordsvtk = LegacyVTKReader(FileNames=['your_directory_path/your_station_coords.vtk'])


    # interpolate to obtain GFs at GPS station and InSAR sampling points
    resampleWithDataset1 = ResampleWithDataset(Input=calculator1, Source=nepal_station_coordsvtk)

    # save data
    SaveData('your_directory_path/%s_topo.csv' % model, proxy=resampleWithDataset1)

    #reset session
    Disconnect()
    Connect()
