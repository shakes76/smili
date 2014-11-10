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
    Joints Scripting Module Global Functions.
    \author Shekhar Chandra

    This module is designed to script batch jobs for automatic segmentation of joints.

    Copyright: (c) 2009 CSIRO, Australia.

    This software is protected by international copyright laws.
    Any unauthorised copying, distribution or reverse engineering is prohibited.

    Licence:
    All rights in this Software are reserved to CSIRO. You are only permitted
    to have this Software in your possession and to make use of it if you have
    agreed to a Software License with CSIRO.

    BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
'''
import os #operating system functions
import fnmatch #wildcard filenaming
import re #reg expressions

'''
Returns a list of integers as strings from filename
'''
def getCaseIDList(filename):
    case_id_strs = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", filename) #extract case number
    return case_id_strs # return list of numbers found as strings
    
'''
Strips the leading zeros and trailing periods from string
'''
def stripID(caseString):
    newCaseString = caseString.rstrip(".")
    if int(newCaseString) != 0:
        newCaseString = newCaseString.lstrip("0")
    return newCaseString
    
'''
Returns the case ID within filename expected at indexed position in name
Case_Index is the index at which the case is expected from a list of numbers in the filename.
Such as the second number in filename etc. This is to cope with multiple numbers in the filename.
'''
def getCaseID(filename, index):
    case_strs = getCaseIDList(filename)
    #~ print case_strs, len(case_strs)
    return int( stripID(case_strs[index]) )
    
'''
Returns the number (as float) within stringValue, removes any symbols and other characters inbetween.
'''
def getNumber(stringValue):
    case_strs = getCaseIDList(stringValue)
    #~ print case_strs, len(case_strs)
    return float( "".join(case_strs) )
    
'''
\brief Return a list of filenames matching the pattern provided.
Pattern is shell like supporting wildcards, such as '*.txt'
'''
def getFileList(path, pattern, prependFullPath=False):
    fileList = [] # init file list
    dirs = os.listdir(path) #pull all filenames from path as strings
    for file in dirs:
        if fnmatch.fnmatch(file, pattern):
            if prependFullPath:
                fileList.append( os.path.abspath(path)+os.sep+file ) 
            else:
                fileList.append(file) 
    return fileList
    
'''
\brief Return a list of filenames matching the pattern provided and a list of cases.
Pattern is shell like supporting wildcards, such as '*.txt'
Case_Index is the index at which the case is expected from a list of numbers in the filename.
Such as the second number in filename etc. This is to cope with multiple numbers in the filename.
'''
def getFileListAndCases(path, case_index, pattern, prependFullPath=False):
    fileList = [] # init file list
    caseList = []
    dirs = os.listdir(path) #pull all filenames from path as strings
    for file in dirs:
        if fnmatch.fnmatch(file, pattern):
            case = getCaseID(file, case_index)
            #~ print case
            caseList.append(case)
            if prependFullPath:
                fileList.append( os.path.abspath(path)+os.sep+file ) 
            else:
                fileList.append(file) 
    return fileList, caseList
    
'''
\brief Return a list of filenames matching the pattern provided and a list of cases.
Cases within ignore list is not within the lists.
Pattern is shell like supporting wildcards, such as '*.txt'
Case_Index is the index at which the case is expected from a list of numbers in the filename.
Such as the second number in filename etc. This is to cope with multiple numbers in the filename.
'''
def getFileListAndCasesWithIgnore(path, ignoreList, case_index, pattern, prependFullPath=False):
    fileList = [] # init file list
    caseList = []
    dirs = os.listdir(path) #pull all filenames from path as strings
    for file in dirs:
        if fnmatch.fnmatch(file, pattern):
            case = getCaseID(file, case_index)
            if case in ignoreList:
                continue
            caseList.append(case)
            if prependFullPath:
                fileList.append( os.path.abspath(path)+os.sep+file ) 
            else:
                fileList.append(file) 
    return fileList, caseList
    
'''
\brief Return a list of filenames matching the pattern and case list provided.
Pattern is shell like supporting wildcards, such as '*.txt'
'''
def getFileListGivenCases(path, cases, pattern, prependFullPath=False):
    fileList = [] # init file list
    dirs = os.listdir(path) #pull all filenames from path as strings
    for file in dirs:
        if fnmatch.fnmatch(file, pattern):
            #~ print file
            case = getCaseID(file, 0)
            if not case in cases:
                continue
            if prependFullPath:
                fileList.append( os.path.abspath(path)+os.sep+file ) 
            else:
                fileList.append(file) 
    return fileList
    
'''
\brief Return a sorted list of filenames matching the pattern provided.
Pattern is shell like supporting wildcards, such as '*.txt'
'''
def getSortedFileList(path, pattern, prependFullPath=False):
    fileList = getFileList(path, pattern, prependFullPath)
    fileList.sort()
    return fileList
    
'''
\brief Return a sorted list of filenames matching the pattern provided.
Pattern is shell like supporting wildcards, such as '*.txt'
Case_Index is the index at which the case is expected from a list of numbers in the filename.
Such as the second number in filename etc. This is to cope with multiple numbers in the filename.
'''
def getSortedFileListAndCases(path, case_index, pattern, prependFullPath=False):
    fileList, caseList = getFileListAndCases(path, case_index, pattern, prependFullPath)
    fileListZipped = zip(fileList, caseList) # zip for easy sorting
    fileListZipped.sort()
    return zip(*fileListZipped) #unzip and return
    
'''
\brief Return a sorted list of filenames matching the pattern provided.
Cases within ignore list is not within the lists.
Pattern is shell like supporting wildcards, such as '*.txt'
Case_Index is the index at which the case is expected from a list of numbers in the filename.
Such as the second number in filename etc. This is to cope with multiple numbers in the filename.
'''
def getSortedFileListAndCasesWithIgnore(path, ignoreList, case_index, pattern, prependFullPath=False):
    fileList, caseList = getFileListAndCasesWithIgnore(path, ignoreList, case_index, pattern, prependFullPath)
    fileListZipped = zip(fileList, caseList) # zip for easy sorting
    fileListZipped.sort()
    return zip(*fileListZipped) #unzip and return
    
'''
\brief Return a sorted list of filenames matching the pattern and case list provided.
Pattern is shell like supporting wildcards, such as '*.txt'
'''
def getSortedFileListGivenCases(path, cases, pattern, prependFullPath=False):
    fileList = getFileListGivenCases(path, cases, pattern, prependFullPath)
    fileList.sort()
    return fileList
    
'''
\brief Write a MilxView Batch (MVB) file for images with case IDs having the name given
'''
def writeImageMVBFile(mvbFileName, nameList, caseList):
    outFile = open(mvbFileName, 'w') #open file for writing

    outFile.write("MILXVIEW_BATCH_FILE\n")
    for name, index in zip(nameList, caseList): # traverse both lists
        outFile.write("CASE_IMAGE " + str(index) + " " + name + "\n")
        
'''
\brief Read list of the ith entries from MilxView Batch (MVB) file
'''
def parseListFile(mapFn, i):
    list = []
    f = open(mapFn,'r')
    lines = f.readlines()
    #Print Lines
    if len(lines) > 0 and 'MILXVIEW_BATCH_FILE' in lines[0]:
        #~ print "Removing MVB identifier"
        lines.pop(0)
    for line in lines:
        #Split the line
        line = line.split()
        #~ print line
        if len(line) > 0:
            #~ print line[i]
            list.append(line[i])
    return list
    
def replace_words(text, word_dic):
    """
    take a text and replace words that match a key in a dictionary with
    the associated value, return the changed text
    """
    rc = re.compile('|'.join(map(re.escape, word_dic)))
    def translate(match):
        return word_dic[match.group(0)]
    return rc.sub(translate, text)
    
def replace_in_file(src, dest, key, newkey):
    '''
    read and replace key with newkey.
    Re-Writing each line instead of copying is a little slow but ok for source files.
    '''
    sourceFile = open(src, "r")
    targetFile = open(dest, "w")
    for line in sourceFile: #for each line
        newline = line
        #handle relative path includes
        if key in line:
            newline = line.replace(key, newkey)
            print "Replace: ", line, " -> ", newline
            
        targetFile.write(newline)
    sourceFile.close()
    targetFile.close()
