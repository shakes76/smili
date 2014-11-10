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
This script concatenates surfaces in case directories
Usage:
/usr/bin/../local/sMILX.0.99/bin/milxModelApp  {--duplicate|-c|--cat
                                        |--colourcat|--split|-s <Scale>|-d
                                        <Decimate>|--smooth <Smooth>
                                        |--laplacian <Laplacian>
                                        |--scalardiff|--scalarstats
                                        |--scalarremove|--scalarcopy|--mse
                                        |--procrustes|--icp
                                        |--thresholdabove <Above>
                                        |--thresholdbelow <Below>|-f <0|1
                                        |2>} [--savetransforms]
                                        [--partitionList] [--rigid]
                                        [--component <Component>] ... 
                                        [--outputformat <Output format>]
                                        [-p <Output Prefix>] [-o <Output>]
                                        [--] [--version] [-h] <Surfaces>
                                        ...


Where: 

   --duplicate
     (OR required)  Simply open and save the input file(s). Work with
     single (-o) or multiple (-p) input
         -- OR --
   -c,  --convert
     (OR required)  Convert the image from current format to the one given
     at output.
         -- OR --
   --cat
     (OR required)  Concatenate surfaces provided
         -- OR --
   --colourcat
     (OR required)  Concatenate surfaces provided and colour them
         -- OR --
   --split
     (OR required)  Split each surface given components
         -- OR --
   -s <Scale>,  --scale <Scale>
     (OR required)  Scale the coordinates of the point
         -- OR --
   -d <Decimate>,  --decimate <Decimate>
     (OR required)  Decimate all the meshes provided using the Quadric
     Decimate algorithm
         -- OR --
   --smooth <Smooth>
     (OR required)  Smooth all the meshes provided using the Windowed Sinc
     Algorithm.
         -- OR --
   --laplacian <Laplacian>
     (OR required)  Smooth all the meshes provided using the Laplacian
     Algorithm.
         -- OR --
   --scalardiff
     (OR required)  Difference in Scalars
         -- OR --
   --scalarstats
     (OR required)  Statistics of Scalars
         -- OR --
   --scalarremove
     (OR required)  Remove the Scalars
         -- OR --
   --scalarcopy
     (OR required)  Copy the Scalars from first mesh to all others while
     removing existing ones.
         -- OR --
   --mse
     (OR required)  Mean Squared Error of Points in models
         -- OR --
   --procrustes
     (OR required)  Similarity alignment of surfaces assuming points have
     correspondence. Use --rigid if rigid alignment is required.
         -- OR --
   --icp
     (OR required)  Iterative Closest Points alignment of surfaces assuming
     points don't have correspondence. Last surface in list is used as the
     reference 'fixed' surface.
         -- OR --
   --thresholdabove <Above>
     (OR required)  Thresold scalars above value
         -- OR --
   --thresholdbelow <Below>
     (OR required)  Thresold scalars below value
         -- OR --
   -f <0|1|2>,  --flip <0|1|2>
     (OR required)  Flip the meshes in the axis provided (0: x-axis, 1:
     y-axis, 2: z-axis).


   --savetransforms
     Save the transformation matrix after surface alignment. For use with
     --icp

   --partitionList
     Partition the list of surfaces provided into two for operating on one
     vs the other

   --rigid
     Rigid alignment of surfaces assuming points have correspondence.

   --component <Component>  (accepted multiple times)
     Surface is a component of the surfaces

   --outputformat <Output format>
     Specify the output format (vtk, vtp, ply, stl, default: same as input)

   -p <Output Prefix>,  --prefix <Output Prefix>
     Output prefix for multiple output

   -o <Output>,  --output <Output>
     Output Surface

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <Surfaces>  (accepted multiple times)
     (required)  Surfaces to operate on


   A diagnostic tool for models/surface operations

Example: 
'''
import filenames
import time
import batch

totalThreads = 6

app = "milxModelApp"

#Constants and  paths
parent_dir = ""
image_path = "images/"
path = "segmentations_all/"
output_path = path
parent_path = filenames.os.getcwd()

#Get the case list from a directory, here the image spaces where the surfaces reside
manualList, caseList = filenames.getSortedFileListAndCases(image_path, 0, '*.hdr', True)
print manualList
print caseList

#create job list to run
jobsList = []
prefix_dirs = []
for file, case in zip(manualList, caseList):
    #output filenames
    prefix_dir = "case_" + str(case) + "/"
    
    #check if results present
    full_prefix_dir = path + prefix_dir + "/"
    if not filenames.os.access(full_prefix_dir, filenames.os.F_OK): #exists? 
        print "No result for", str(case)
        continue
    objList = filenames.getSortedFileList(full_prefix_dir, '*_Relaxed_Combined.nii.gz', True)
    if not objList:
        print "Not result present at", full_prefix_dir
        continue
    objList1 = filenames.getSortedFileList(full_prefix_dir, '*_Relaxed_object_001.vtk', True)
    if not objList1:
        print "Not result present at", full_prefix_dir
        continue
    objList2 = filenames.getSortedFileList(full_prefix_dir, '*_Relaxed_object_002.vtk', True)
    if not objList2:
        print "Not result present at", full_prefix_dir
        continue
    surface1 = objList1[0]
    surface2 = objList2[0]
    
    name, ext = filenames.os.path.basename(objList[0]).split('.', 1)
    output_name = output_path + prefix_dir + name + ".vtk"
    
    command = app + " --cat " + surface1 + " " + surface2 + " -o " + output_name
    print >> batch.stderr, command # Print to std err
    
    commandArgs = []
    commandArgs.append("--cat")
    commandArgs.append(surface1)
    commandArgs.append(surface2)
    commandArgs.append("-o")
    commandArgs.append(output_name)
    #~ commandArgs.append("--padding")
    #~ commandArgs.append("0")
    job = {
        'application' : app,
        'arguments' : commandArgs,
        'identifier' : "Case " + str(case)
        }
    jobsList.append(job)

#~ print jobsList

#For each sorted filename, compute the segmentation of the image
start = time.time()
#run the jobs over multiple threads
batch.processJobs(jobsList, totalThreads, True)
end = time.time()
elapsed = end - start
print "Transforms took " + str(elapsed) + " secs or " + str(elapsed/60) + " mins in total"
