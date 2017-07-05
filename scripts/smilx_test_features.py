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
This script setups a scene from multiple files and outputs a screenshot.
The scene consists of a vector field and an image slice with a predefined camera position.
The vector field is generated first from a deformation field vector image.
To use it, run smilx, in python console type "execfile('smilx_screenshot.py')"
'''
#names of images, set them to empty string to get popup dialogs
path = 'E:/Dev/smili-release/scripts/'
imageName = path+"data/Atlas_MRI_Mean2_R_Preprocessed_x2.nii.gz"
labelName = path+"data/Atlas_MRI_Mean2_R_labels_x2.nii.gz"
surfaceName1 = path+"data/Atlas_MRI_surface_R_fem.vtk"
surfaceName2 = path+"data/Atlas_MRI_surface_R_ace.vtk"

MainWindow.loadFile(imageName) # load image
image1 = MainWindow.activeImage() # get pointer to image window

MainWindow.loadFile(imageName) # load image
image2 = MainWindow.activeImage() # get pointer to image window

MainWindow.loadFile(labelName) # load image
labelImage = MainWindow.activeImage() # get pointer to image window

MainWindow.loadFile(surfaceName1) # load image
surface1 = MainWindow.activeModel() # get pointer to image window
MainWindow.loadFile(surfaceName2) # load image
surface2 = MainWindow.activeModel() # get pointer to image window

surface1.importViewFrom(image1) #import the slice view to the vector field
#~ surface2.importViewFrom(surface1) #import the slice view to the surface window

image1.gradientMagnitude() #edges
image2.overlay(labelName)

MainWindow.tileTab()

#~ MainWindow.link() #link all the windows
MainWindow.saveScreen(path+'screenie.png') #snapshot
