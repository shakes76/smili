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
Match the info of one surface with another. Save result to file.
'''
import os, re
 
meanFile = "mean.vtk"
meanPath = "/home/cha588/dgv-py-test/"
 
meshPath = "/home/cha588/data/ave_to_MRs_proc_meshes_nrr/"
 
outputPath = "/home/cha588/dgv-py-test/init_mean_meshes/"
outputPrefix = "init_mean_"
outputExt = ".vtk"
 
dirs = os.listdir(meshPath)
for file in dirs:
    MainWindow.loadFile(meanPath + meanFile) # load mean mesh
 
    mdl_mean.matchInfo(meshPath + file) #Match the mean to other mesh
     
    #find case id in filename and remember it
    case_id_str = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", file)
    case = int( case_id_str[0] )
     
    outputName = "init_mean_" + str(case) + outputExt
 
    milxQtFile.saveModel(outputPath+outputName, mdl_mean) # save result
 
    mdl_mean.close() # close (since its delete on close
     
print "Matching Done"
