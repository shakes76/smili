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
Batch job processing module.

A joblist is passed into the processJobs() function, and it processes the jobs ensuring all threads are always busy. Number of Threads native to the system is available as the variable 'cores'.

The job list consists of a list or a list of dictionaries:
A job list simply has a list command strings which are executed in the order they appear.
The job list dictionary has to have the keys 'application', 'identifier' and 'arguments' defined (see usage example at end).
Identifier key is used for output purposes and is not essential. It can be left blank.
Arguments needs to be list of strings, not one long string with spaces.

Notes:
home '~' is not recognised by the Python subprocess module. In the latter case, qualify the full path.
'''
from Queue import Queue
from threading import Thread
from sys import stderr
import subprocess
from multiprocessing import cpu_count #for cpu counting

cores = cpu_count()

#----------------------------------------------------------------------------

def processJobs(jobs, numberOfThreads, verboseMode=False):
    ''' 
    Executes the list of jobs as n queued jobs
    jobs - a list with dictionary entries of job(s) to run. It is expected to be a dictionary with 'application', 'identifier' and 'arguments' keys
    Arguments needs to be list of strings, not one long string with spaces.
    numberOfThreads - how many threads to use
    
    USAGE:
    import batch
        
    jobsToRun = []

    job1 = {
        'application' : "echo",
        'arguments' : ["'Yay'"], #must be a list
        'identifier' : "Yay Job"}
    job2 = {
        'application' : "echo",
        'arguments' : ["'Woo hoo'"], #must be a list
        'identifier' : "Woohoo Job"}
    job3 = {
        'application' : "echo",
        'arguments' : ["'Yippee'"], #must be a list
        'identifier' : "Yippee Job"}
        
    jobsToRun.append(job1)
    jobsToRun.append(job2)
    jobsToRun.append(job3)
    batch.processJobs(jobsToRun, batch.cores)
    '''
    
    def dictWorker():
        ''' 
        A standard worker thread - Each iteration pulls a job out of the queue and runs it.
        '''
        while True:
            item = jobQueueDict.get()
            # So spawn processes here and wait for them to finish
            application = item['application']
            arguments = item['arguments']
            identifier = item['identifier']
            
            #create command list
            command = []
            command.append(application.strip()) #strip used to ensure no spaces
            for arg in arguments:
                command.append(arg)
            
            if verboseMode:
                print >> stderr, "Started: ", identifier # Print to std err
            subprocess.call(command) # run command
            if verboseMode:
                print >> stderr, "Done: ", identifier # Print to std err
            jobQueueDict.task_done()

    jobQueueDict = Queue()
    for i in range(numberOfThreads):
         t = Thread(target=dictWorker)
         t.daemon = True
         t.start()
    
    for item in jobs:
        jobQueueDict.put(item)

    jobQueueDict.join()       # block until all tasks are done

#----------------------------------------------------------------------------

def processJobsFromList(jobs, numberOfThreads, verboseMode=False):
    ''' 
    Executes the list of jobs as n queued jobs
    jobs - a list of job(s) to run. Each entry is expected to be the command string to execute
    numberOfThreads - how many threads to use
    
    USAGE:
    import batch
        
    jobsToRun = []

    job1 = "echo 'Yay'"
    job2 = "echo 'Woo hoo'"
    job3 = "echo 'Yippee'"
        
    jobsToRun.append(job1)
    jobsToRun.append(job2)
    jobsToRun.append(job3)
    batch.processJobsFromList(jobsToRun, batch.cores)
    '''
    
    def listWorker():
        ''' 
        A standard worker thread - Each iteration pulls a job out of the queue and runs it.
        '''
        while True:
            item = jobQueueList.get()
            # So spawn processes here and wait for them to finish
            tokens = item.split()
            
            #create command list
            command = []
            for arg in tokens:
                command.append(arg.strip())  #strip used to ensure no spaces
            
            if verboseMode:
                print >> stderr, "Start" # Print to std err
            subprocess.call(command) # run command
            if verboseMode:
                print >> stderr, "Done" # Print to std err
            jobQueueList.task_done()

    jobQueueList = Queue()
    for i in range(numberOfThreads):
         t = Thread(target=listWorker)
         t.daemon = True
         t.start()
    
    for item in jobs:
        jobQueueList.put(item)

    jobQueueList.join()       # block until all tasks are done

def processJobsFromListInShell(jobs, numberOfThreads, verboseMode=False):
    ''' 
    Executes the list of jobs as n queued jobs
    jobs - a list of job(s) to run. Each entry is expected to be the command string to execute
    numberOfThreads - how many threads to use
    
    USAGE:
    import batch
        
    jobsToRun = []

    job1 = "echo 'Yay'"
    job2 = "echo 'Woo hoo'"
    job3 = "echo 'Yippee'"
        
    jobsToRun.append(job1)
    jobsToRun.append(job2)
    jobsToRun.append(job3)
    batch.processJobsFromList(jobsToRun, batch.cores)
    '''
    
    def listWorker():
        ''' 
        A standard worker thread - Each iteration pulls a job out of the queue and runs it.
        '''
        while True:
            item = jobQueueList.get()
            # So spawn processes here and wait for them to finish
            tokens = item.split()
            
            #create command list
            command = []
            for arg in tokens:
                command.append(arg.strip())  #strip used to ensure no spaces
            
            if verboseMode:
                print >> stderr, "Start" # Print to std err
            subprocess.Popen(item, shell=True).wait # run command
            if verboseMode:
                print >> stderr, "Done" # Print to std err
            jobQueueList.task_done()

    jobQueueList = Queue()
    for i in range(numberOfThreads):
         t = Thread(target=listWorker)
         t.daemon = True
         t.start()
    
    for item in jobs:
        jobQueueList.put(item)

    jobQueueList.join()       # block until all tasks are done
