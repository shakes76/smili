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
This script opens a pv map and its corresponding image, thresholds and overlays the two, loads camera view, links windows and disables interpolation.
To use it, run smilx, in python console type "execfile('smilx_pv_visualisation.py')"
'''
#names of images, set them to empty string to get popup dialogs
imageName = "PAVEL_FLAWS_INV2.nii.gz"
pvImageName = "PAVEL_BiExp_Combined_PV_GM.nii.gz"
tmpOutName = "thresholded_pvimage.nii.gz"

MainWindow.loadFile(imageName) # load image
image = MainWindow.activeImage() # get pointer to image window

MainWindow.loadFile(pvImageName) # load image
pvImage = MainWindow.activeImage() # get pointer to image window

#process pv map for display
belowLevel = 0.01
aboveLevel = 0.25
pvImage.binaryThreshold(255, belowLevel, aboveLevel)

milxQtFile.saveImage(tmpOutName, pvImage) # save result

#overlay the pv processed result
pvImage.viewToSagittal()
pvImage.loadView("camera.cam")
pvImage.interpolateDisplay()
image.overlay(tmpOutName)
image.viewToSagittal()
image.loadView("camera.cam")
image.interpolateDisplay()

MainWindow.tileTab()
MainWindow.link() #link all the windows
