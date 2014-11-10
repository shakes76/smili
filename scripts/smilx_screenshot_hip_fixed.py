#!/usr/bin/smilx
'''=========================================================================
  The Software is copyright (c) Commonwealth Scientific and Industrial Research Organisation (CSIRO)
  ABN 41 687 119 230.
  All rights reserved.

  Licensed under the CSIRO BSD 3-Clause License
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License in the file LICENSE.md or at

  https://stash.csiro.au/projects/SMILI/repos/smili/browse/license.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
========================================================================='''
'''
This script generates screenshot from a set of whole meshes and clipped meshes (in the same space).
Each model is loaded and put into the same view and overlaid. The screenshot is then generated.
The meshes in view are fixed only scalars change
'''
import os, re

execfile('filenames.py')
 
mesh_path = "./"
clip_path = "clipped/"

objects = ['fem', 'ace']
maxValues = [ [1.75, 5.0, 1.0], [1.75, 7.5, 5.0]]
cameras = ['fem.cam', 'ace.cam'] #label values
#~ methods = [ 'robust', 'robustregions', 'weight', 'weightregions', 'weightregionsfast', 'std', 'dynamic' ]
methods = [ 'robust', 'weight', 'weightregions', 'weightregionsfast' ]
stats = [ 'Mean per Vertex', 'Maximum per Vertex', 'Variance per Vertex' ]
stats_title = [ 'Mean Error (mm)', 'Maximum Error (mm)', 'Error Variance (mm)' ]
stats_name = [ '_Mean', '_Maximum', '_Variance' ]
 
outputPath = "figures/"
outputExt = ".png"
outputPrefix = "fig_"

MainWindow.working(-1) #for progress bar
for obj, max, cam in zip(objects, maxValues, cameras):
    meshList = getSortedFileList(mesh_path, obj+'_*.vtp', True)
    clipList = getSortedFileList(clip_path, 'clip_'+obj+'_*.vtp', True)
    print meshList
    print clipList

    for mesh, clip in zip(meshList, clipList):
        MainWindow.loadFile(meshList[0]) # load model
        currentModelWhole = MainWindow.activeModel() # get the loaded model
        currentModelWhole.background(True) #white
        currentModelWhole.orientDisplay() #disable human figure
        currentModelWhole.enableInterpolateDisplay() 
        currentModelWhole.removeScalars() 
        
        MainWindow.loadFile(clipList[0]) # load model
        currentModel = MainWindow.activeModel() # get the surface
        currentModel.loadScalars(clip) # load scalars from different mesh, possible because have correspondence
        currentModel.background(True) #white
        currentModel.orientDisplay() #disable human figure
        currentModel.enableInterpolateDisplay() # phong shading
        currentModel.colourMapToHOT()
        
        MainWindow.tileTab() #tile cause camera view saved is optimised for it
        
        currentModel.loadView(cam)
        currentModelWhole.loadView(cam)
        
        for stat, title, stat_name, maxVal in zip(stats, stats_title, stats_name, max):
            currentModel.disableScale() #stop it from asking auto scalarbar
            currentModel.showArray(stat)
            currentModel.enableScale(title, True, 0, maxVal)
            currentModelWhole.importFrom(currentModel)
         
            outputName = outputPath + outputPrefix + os.path.splitext( os.path.basename(mesh) )[0] + stat_name + outputExt
            print "Saving ", outputName
            MainWindow.setActiveWindow(currentModelWhole) #ensure right window is active to save screenshot from
            MainWindow.saveScreen(outputName) # save result
         
        currentModel.close() # close (since its delete on close)
        currentModelWhole.close() # close (since its delete on close)
        
        #~ break #for testing
    
    #~ break #for testing
MainWindow.done(-1)
     
print "Processing Done"
