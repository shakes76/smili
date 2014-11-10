#!/usr/bin/python
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
This script generates 3D rendering screen shots of the results
Here the app is a overlay app with the syntax:
Usage:
    milxOverlay  {--axial|--coronal
                                        |--saggital} [--nohuman]
                                        [--wireframes] [--wireframe]
                                        [--autocolour] [--noslices]
                                        [--outline] [--removescalars]
                                        [--nosaggital] [--nocoronal]
                                        [--noaxial] [--inverse] [--white]
                                        [--onscreen] [--loadviewfile <Load
                                        View File>] [--loadview] [--HOT]
                                        [--AAL] [--NIH_FIRE] [--NIH]
                                        [--rainbow] [--vtk] [--jet]
                                        [--saggitalslice <Saggital Slice>]
                                        [--coronalslice <Coronal Slice>]
                                        [--axialslice <Axial Slice>] [-x
                                        <Width>] [-y <Height>]
                                        [--setscalars <Set Scalars>]
                                        [--loadscalars <Scalars>]
                                        [--isovalue <Isovalue>] [--opacity
                                        <Opacity>] [--bluebackground <Blue
                                        Component>] [--greenbackground
                                        <Green Component>] [--redbackground
                                        <Red Component>] [--max <Maximum
                                        Scalar>] [--min <Minimum Scalar>]
                                        [--isosurface <Isosurface>] [-t
                                        <Transform>] [-l <Label>] [-s
                                        <Scalar Mask>] [-v <Vectors>] [-i
                                        <Image>] [-o <Output>] [--]
                                        [--version] [-h] <Surfaces> ...

Where: 

   --axial
     (OR required)  Show axial view based on position of first surface.
         -- OR --
   --coronal
     (OR required)  Show coronal view based on position of first surface.
         -- OR --
   --saggital
     (OR required)  Show saggital view based on position of first surface.


   --nohuman
     Disable human orientation glyph.

   --wireframes
     Display all surfaces as a wireframe models.

   --wireframe
     Display initial surface as a wireframe model.

   --autocolour
     Colour each additional model automatically based on order.

   --noslices
     Prevent any slices from being shown. Other elements from image option
     etc. will still be displayed.

   --outline
     Display outline box for the model. Model outline overided by image if
     present.

   --removescalars
     Remove the scalars of the models.

   --nosaggital
     Do not show saggital slice when image is given.

   --nocoronal
     Do not show coronal slice when image is given.

   --noaxial
     Do not show axial slice when image is given.

   --inverse
     Use the inverse transform when transform file is provided.

   --white
     Make background white rather than default gradient colour.

   --onscreen
     Enable on screen rendering, i.e. display the rendering as an
     interactive window.

   --loadviewfile <Load View File>
     Load saved view from file (use onscreen mode to save view files)

   --loadview
     Load saved view (use smilx or onscreen render mode to view and save
     with Right Click->View->Save View

   --HOT
     Change colourmap of the scalars to HOT

   --AAL
     Change colourmap of the scalars to AAL

   --NIH_FIRE
     Change colourmap of the scalars to NIH Fire

   --NIH
     Change colourmap of the scalars to NIH

   --rainbow
     Change colourmap of the scalars to the rainbow map

   --vtk
     Change colourmap of the scalars to blue-red (rainbow VTK) map

   --jet
     Change colourmap of the scalars to the Jet map

   --saggitalslice <Saggital Slice>
     Set the saggital slice of the image.

   --coronalslice <Coronal Slice>
     Set the coronal slice of the image.

   --axialslice <Axial Slice>
     Set the axial slice of the image.

   -x <Width>,  --width <Width>
     Set the width of the window.

   -y <Height>,  --height <Height>
     Set the height of the window.

   --setscalars <Set Scalars>
     Set active scalars of model as array name given

   --loadscalars <Scalars>
     Load scalars for the mesh from another mesh

   --isovalue <Isovalue>
     Set the label or iso-surface value for option --isosurface.

   --opacity <Opacity>
     Set the opacity of the models

   --bluebackground <Blue Component>
     Set the blueness value (0-1) for the background.

   --greenbackground <Green Component>
     Set the greenness value (0-1) for the background.

   --redbackground <Red Component>
     Set the redness value (0-1) for the background.

   --max <Maximum Scalar>
     Set the maximum value for scalars on model.

   --min <Minimum Scalar>
     Set the minimum value for scalars on model.

   --isosurface <Isosurface>
     Generate Iso-surface (surface from label image) from the image

   -t <Transform>,  --transform <Transform>
     Transform (ITK Format) to apply to objects being rendered

   -l <Label>,  --label <Label>
     Overlay binary image/label on the image

   -s <Scalar Mask>,  --scalarmask <Scalar Mask>
     Clip the main mesh based on binary scalar values in model

   -v <Vectors>,  --vectors <Vectors>
     Generate and display vector field from model. If no vectors are
     present then scalars and normals will be used

   -i <Image>,  --image <Image>
     Load and generate image slices as planes (axial, coronal & saggital)

   -o <Output>,  --output <Output>
     Output Screenshot name

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <Surfaces>  (accepted multiple times)
     (required)  Surfaces to overlay

   An overlay tool for models/images

Example: milxOverlay /home/cha588/data/Hip_new/workspace/pseudo_unilaterals/segmentations/case_1//hip_1_001_Segmentation_object_001_Clipped.vtk /home/cha588/data/Hip_new/workspace/pseudo_unilaterals/segmentations/case_1//hip_1_001_Segmentation_object_002_Clipped.vtk -o /home/cha588/data/Hip_new/workspace/pseudo_unilaterals/segmentations/hip_pseudo_unilaterals_1.png --image /home/cha588/data/Hip_new/workspace/pseudo_unilaterals/segmentations/case_1//hip_1__Propagated_Image_to_Atlas_000.nii.gz --white --loadview --nocoronal --noaxial --wireframes --width 600 --height 800

'''
import filenames
import time
import batch

#~ totalThreads = batch.cores/2
totalThreads = 2

app = "milxOverlay"

options = "--white --wireframe --coronal --autocolour "
windowsize = "--width 640 --height 480 "

#Constants and  paths
object = "rectum"
mesh_path1 = 'output_rectum_vlow/'
mesh_path2 = mesh_path1
output_path = mesh_path1
output_prefix = 'screen_initial_' + object + '_'

#The ordering of the surfaces in the SSM are assumed to be alphabetical
meshList1, caseList = filenames.getSortedFileListAndCases(mesh_path1, 0, '*_RawSurfaces_*.vtp', True)
meshList2 = filenames.getSortedFileListGivenCases(mesh_path2, caseList, '*_InitialSurfaces_*.vtp', True)
print meshList1
print meshList2
print caseList

#command specific options
other_options = " "

command_common_options = options + windowsize + other_options
commonList = command_common_options.split() #create indexable list of above string

print "Command Options: "
print commonList

#create job list to run
jobsList = []
for mesh1, mesh2, case in zip(meshList1, meshList2, caseList):
    
    #Overlay output
    output_name = output_path + output_prefix + str(case).zfill(3) + ".png"
    
    #case command options
    prefix_option = "-o " + output_name + " "
    command_options = prefix_option
    commandList = command_options.split() #create indexable list of above string
    
    command = app + " " + mesh1 + " " + mesh2 + " " + command_options + command_common_options
    print >> batch.stderr, command # Print to std err
    
    commandArgs = []
    commandArgs.append(mesh1)
    commandArgs.append(mesh2)
    for option in commandList:
        commandArgs.append(option)
    for option in commonList:
        commandArgs.append(option)
    job = {
        'application' : app,
        'arguments' : commandArgs,
        'identifier' : "Case " + str(case)
        }
    jobsList.append(job)
    #~ break

#~ print jobsList

#For each sorted filename, compute the segmentation of the image
start = time.time()
#run the jobs over multiple threads
batch.processJobs(jobsList, totalThreads, True)
end = time.time()
elapsed = end - start
print "Overlays took " + str(elapsed) + " secs or " + str(elapsed/60) + " mins in total"
