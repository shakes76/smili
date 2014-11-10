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
imageName = "RegOut_7.nii.gz"
defImageName = "RegDFOut_7.nii.gz"

MainWindow.loadFile(imageName) # load image
image = MainWindow.activeImage() # get pointer to image window

MainWindow.loadFile(defImageName) # load image
defImage = MainWindow.activeImage() # get pointer to image window

defImage.vectorField(8, 3.0) #subsampling and scaling given

MainWindow.tileTab()

defField = MainWindow.activeModel() # get pointer to image window
defField.loadView("camera.cam")
defField.importViewFrom(image) #import the slice view to the vector field

#~ MainWindow.link() #link all the windows
MainWindow.saveScreen("screenie.png") #snapshot
