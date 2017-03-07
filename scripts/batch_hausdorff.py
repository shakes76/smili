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
This script generates Hausdorff distances for the hip results
Here the app is a distance app with the syntax:
Usage:
   milxHausdorffDistance  [-c <Case>] [-s] --labelvalue <Label Value>
                                -l <Label> [-p <Output Prefix>] [-o
                                <Output>] [--] [--version] [-h] <Surfaces>
                                ...

Where: 

   -c <Case>,  --case <Case>
     Set the case ID being done. Used to name extra output.

   -s,  --symmetric
     Compute forward and backward distances. This is required to get
     Hausdorff distance.

   --labelvalue <Label Value>
     (required)  Set the label value for option --label.

   -l <Label>,  --label <Label>
     (required)  Compute the distances from the labelled image to the
     surface(s) provided.

   -p <Output Prefix>,  --prefix <Output Prefix>
     Output prefix for multiple output

   -o <Output>,  --output <Output>
     Output model name

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <Surfaces>  (accepted multiple times)
     (required)  Surfaces to compute the distances with.


   A Hausdorff Distance tool for models

Example: 

'''
import filenames
import time
import batch

totalThreads = batch.cores/2

#Constants and  paths
#~ parent_dir = filenames.os.getcwd()+'/'
parent_dir = ''
manual_path = 'manuals_renamed/'
output_path = 'Hausdorff/'
home_dir = filenames.os.path.expanduser("~")
smili = home_dir+'/Dev/smili/build/'
app = smili+"bin/milxHausdorffDistance"

options = " "
#~ options = " --symmetric "

#~ result_dirs = ['segmentations', 'segmentations_dyn', 'segmentations_robust', 'segmentations_weight', 'segmentations_weight_fast']
#~ result_dirs_names = ['std', 'dyn', 'robust', 'weight', 'weightfast']
result_dirs = ['results_clipped']
result_dirs_names = ['']
objects = ['bladder', 'rectum', 'prostate', 'prostate_T2MR']
object_case_index = [0, 0, 0, 1]
objects_values = ['1', '1', '1', '1']

for dir, dirName in zip(result_dirs, result_dirs_names):
    input_path = parent_dir + dir + '/'
    for object, value, case_index in zip(objects, objects_values, object_case_index):
        object_output_path = output_path+object+'/'
        output_prefix = object + "_" + dirName + "_"
        
        if not filenames.os.access(object_output_path, filenames.os.F_OK): #exists? No then create
            filenames.os.mkdir(object_output_path)
            
        #is manuals present, else dont bother
        manualList = filenames.getSortedFileList(manual_path+object, '*.nii.gz', True)
        if not manualList:
            print "Dataset doesn't have manuals. Skipping."
            continue
            
        #The ordering of the surfaces in the SSM are assumed to be alphabetical
        manualList, manCaseList = filenames.getSortedFileListAndCases(manual_path+object, case_index, '*.nii.gz')
        print manualList
        print manCaseList

        commonList = options.split() #create indexable list of above string

        print "Command Options: "
        print commonList

        #create job list to run
        jobsList = []
        prefix_dirs = []
        outNames = []
        for file, case in zip(manualList, manCaseList):
            #output filenames
            prefix_dir = object
            full_prefix_dir = input_path + prefix_dir+'/'
            print "Case", case, ":", file
            
            #check if result present
            objList, cases = filenames.getSortedFileListAndCases(full_prefix_dir, case_index, 'asm_'+object+'_*' + '.vtk')
            if not case in cases:
                print "Not result for", case,"present at", full_prefix_dir
                continue
            
            manual = manual_path+object+'/'+file
            index = cases.index(case)
            result = full_prefix_dir+objList[index]
            print "Result found:", result
                
            #Hausdorff out name
            output_name = object_output_path + output_prefix + str(case).zfill(3) + ".vtk"
            outNames.append(output_name)
            
            #case command options
            image_option = "--label " + manual + " "
            value_option = "--labelvalue " + str(value) + " "
            out_option = "-o " + output_name + " "
            case_option = "-c " + str(case).zfill(3) + " "
            prefix_option = "-p " + object_output_path + output_prefix + " "
            command_options = prefix_option + out_option + image_option + value_option + case_option
            commandList = command_options.split() #create indexable list of above string
            
            command = app + " " + result + " " + command_options + options
            print command
            jobsList.append(command)
            
            #~ break

        #~ print jobsList

        #For each sorted filename, compute the segmentation of the image
        start = time.time()
        #run the jobs over multiple threads
        batch.processJobsFromList(jobsList, totalThreads, True)
        end = time.time()
        elapsed = end - start
        print "Hausdorffs took " + str(elapsed) + " secs or " + str(elapsed/60) + " mins in total"
