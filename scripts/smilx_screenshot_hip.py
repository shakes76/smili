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
'''
import os, re

execfile('filenames.py')
 
mesh_path = "./"
clip_path = "clipped/"

objects = ['fem', 'ace']
maxValues = [1.5, 2.0]
cameras = ['fem.cam', 'ace.cam'] #label values
#~ methods = [ 'robust', 'robustregions', 'weight', 'weightregions', 'weightregionsfast', 'std', 'dynamic' ]
methods = [ 'robust', 'weight', 'weightregions', 'weightregionsfast' ]
 
outputPath = "figures/"
outputExt = ".png"

MainWindow.working(-1) #for progress bar
for obj, max, cam in zip(objects, maxValues, cameras):
    meshList = getSortedFileList(mesh_path, obj+'_*.vtp', True)
    clipList = getSortedFileList(clip_path, 'clip_'+obj+'_*.vtp', True)
    print meshList
    print clipList

    for mesh, clip in zip(meshList, clipList):
        outputPrefix = "fig_"
        outputName = outputPath + outputPrefix + os.path.splitext( os.path.basename(mesh) )[0] + outputExt
        
        MainWindow.loadFile(mesh) # load image
        currentModelWhole = MainWindow.activeModel() # get the loaded model
        currentModelWhole.background(True)
        currentModelWhole.orientDisplay() #disable human figure
        currentModelWhole.enableInterpolateDisplay() 
        currentModelWhole.removeScalars() 
        
        MainWindow.loadFile(clip) # load image
        currentModel = MainWindow.activeModel() # get the surface
        currentModel.background(True)
        currentModel.orientDisplay() #disable human figure
        currentModel.enableInterpolateDisplay()
        currentModel.colourMapToHOT()
        currentModel.enableScale("Mean Error", True, 0, max)
        
        MainWindow.tileTab()
        
        currentModel.loadView(cam)
        currentModelWhole.loadView(cam)
        
        currentModelWhole.importFrom(currentModel)
         
        print "Saving ", outputName
        MainWindow.setActiveWindow(currentModelWhole) #ensure right window is active to save screenshot from
        MainWindow.saveScreen(outputName) # save result
         
        currentModel.close() # close (since its delete on close)
        currentModelWhole.close() # close (since its delete on close)
        
        #~ break #for testing
    
    #~ break #for testing
MainWindow.done(-1)
     
print "Processing Done"
