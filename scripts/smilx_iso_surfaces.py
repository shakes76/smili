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
This script opens a model one at a time from the path provided and processes it. The processing is
Clean, Smooth, Decimate by 25%, Smooth, Decimate by 25%, Smooth.
The resulting model is then saved to disk with naming convention provided.
'''
import os, re

execfile('filenames.py')
 
manual_path = "combined/"

objects = ['fem', 'ace']
objects_values = [30, 90] #label values
 
outputPath = "meshes/"
outputExt = ".vtk"

manualList, caseList = getSortedFileListAndCases(manual_path, 0, '*.nii.gz', True)
print manualList

MainWindow.working(-1) #for progress bar
for obj, label in zip(objects, objects_values):

    for manual, case in zip(manualList, caseList):
        outputPrefix = "manual_surface_" + obj + "_"
        outputName = outputPath + outputPrefix + str(case) + outputExt
        
        if os.path.exists(outputName):
            print "Found", outputName, "Skipping."
            continue
        
        MainWindow.loadFile(manual) # load image
         
        currentImg = MainWindow.activeImage() # get the loaded image
        currentImg.hide()
        
        #threshold and iso surface
        currentImg.binaryThreshold(1, label, label)
        currentImg.surface(0.5)
        currentImg.hide()
        
        currentModel = MainWindow.activeModel() # get the surface
     
        #process model
        currentModel.clean()
        currentModel.quadricDecimate(0.5)
         
        print "Saving ", outputName
        milxQtFile.saveModel(outputName, currentModel) # save result
         
        currentModel.close() # close (since its delete on close)
        currentImg.close() # close (since its delete on close)
        
        #~ break #for testing
    
    #~ break #for testing
MainWindow.done(-1)
     
print "Processing Done"
