
# This script generate the illustration pictures for the documentation.
#
# To run this script, you have to install paraview, see:
#
# http://www.paraview.org/paraview/resources/software.php
#
# Eventually, set up the paths (figure_path, tutorial_data_path) according to your own installation.
# (The default values should be ok.)
#
# Make sure that pvpython is in your path of executables.
#
# After all the tutorial programs have been executed, run the following
# command in the same directory:
#
#   pvpython generate_doc_figures.py Documentation/Figure
#

from paraview.simple import *
from os import remove, mkdir, curdir
from os.path import join, isdir
from sys import argv, exit

# we need at least the output directory
if len(argv) <= 1:
	exit('Synopsis: pvpython generate_doc_figures.py dest-dir [src-dir]')

figure_path = argv[1]

# default for the input directory is the current one
if len(argv) <= 2:
    tutorial_data_path = curdir
else:
    tutorial_data_path = argv[2]

collected_garbage_file = []

if not isdir(figure_path):
    mkdir(figure_path)
    
## [tutorial5]
# tutorial 5
exts = [".vtu",".dat"]
for case in range(0,903):

    for ext in exts:
        data_file_name = join(tutorial_data_path, "blackoil-output-"+"%(case)d"%{"case": case}+ext)   
        collected_garbage_file.append(data_file_name)
 

cases = ["0", "300", "600", "900"]
for case in cases:
    data_file_name = join(tutorial_data_path, "blackoil-output-"+case+".vtu")
    grid = XMLUnstructuredGridReader(FileName = data_file_name)
    grid.UpdatePipeline()
    Show(grid)
    dp = GetDisplayProperties(grid)
    dp.Representation = 'Surface'
    dp.ColorArrayName = 'sat[1]'
    sat = grid.CellData.GetArray(1)
    sat_lookuptable = GetLookupTableForArray( "sat[1]", 1, RGBPoints=[0, 1, 1, 1, 1, 0, 0, 1])
    dp.LookupTable = sat_lookuptable
    view = GetActiveView()
    view.Background = [1, 1, 1]
    camera = GetActiveCamera()
    camera.SetPosition(0, 4000, 500)
    camera.SetViewUp(0, 0, -1)
    camera.SetViewAngle(30)
    camera.SetFocalPoint(1000, 0, 3000)
    Render()
    WriteImage(join(figure_path, "blackoil-output-"+case+".png"))
Hide(grid)
## [tutorial5]

# remove temporary files
for f in collected_garbage_file:
    remove(f)

