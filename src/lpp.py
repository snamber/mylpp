oneline = "writing pp-data in vtk format automatically, saving memory"

docstr = """this is the docstr of LIGGGHTSPostProcessing"""

from dump import dump
from math import floor
from math import ceil
import vtk
import glob
import multiprocessing
import sys
import time
import os
import exceptions
# changes py p.s. - include a few command line options
import getopt

class lpp:
  
  #=============================================================================
  # creates a filelist, seperates it to sublists
  # creates multiple processes
  # calls lppWorker for each sublist
  #   lppWorker
  #   calls dump, vtk and manyGran for the given list of files
  #   returns 0
  #=============================================================================
  
  def __init__(self, *list, **kwargs):
    # do argument parsing, raise errors if non-integers were given
    # this can be changed if one wants less overhead but use more memory:
    # make this figure higher but if possible a multiple of 8
    self.cpunum      = multiprocessing.cpu_count()
    self.chunksize   = 8
    self.overwrite   = True

    if "--chunksize" in kwargs: 
      try:
        if int(kwargs["--chunksize"]) > 0:
          self.chunksize = int(kwargs["--chunksize"])
        else: raise ValueError
      except ValueError:
        raise ValueError, "Invalid or no argument given for chunksize"

    if "--cpunum" in kwargs: 
      try:
        if int(kwargs["--cpunum"]) > 0 and int(kwargs["--cpunum"]) <= self.cpunum:
          self.cpunum = int(kwargs["--cpunum"])
        else: raise ValueError
      except ValueError:
        raise ValueError, "Invalid or no argument given for cpunum"

    # do not overwrite existing files
    if "--no-overwrite" in kwargs:
      self.overwrite = False
    
    # suppress output with 'False'
    if "--debug" in kwargs: self.debugMode = True
    else: self.debugMode = False
    
    if "--quiet" in kwargs:
      self.output = False
      self.debugMode = False
    else: self.output = True

    if self.output:
      print "starting LIGGGHTS memory optimized parallel post processing"
      print "chunksize:", self.chunksize, "-->",self.chunksize,\
        "files are processed per chunk. If you run out of memory reduce chunksize."
    starttime = time.time()    

    if self.debugMode: print "number of process:", os.getpid()
    
    # check whether file-list is nonempty
    self.flist = list[0]
    listlen = len(self.flist)
    if listlen == 0 and len(list) == 1:
      raise StandardError, "no dump file specified"
    if listlen == 1 and self.overwrite == False:
      raise StandardError, "Cannot process single dump files with --no-overwrite."
    
    if self.output:
      print "Working with", self.cpunum, "processes..."
    
    # seperate list in pieces+rest
    self.slices = []
    
    residualPresent = int(bool(listlen-floor(listlen/self.chunksize)*self.chunksize))
    
    for i in xrange(int(floor(listlen/self.chunksize))+residualPresent):
      slice = self.flist[i*self.chunksize:(i+1)*self.chunksize]
      self.slices.append(slice)
    self.flist = []

    output = ""
    if "-o" in kwargs: output = kwargs["-o"]
    # generate input for lppWorker
    dumpInput = [{"filelist":self.slices[i],\
      "debugMode":self.debugMode,\
      "output":output,\
      "overwrite":self.overwrite} \
      for i in xrange(len(self.slices))]
    
    if self.debugMode: print "dumpInput:",dumpInput
    
    numberOfRuns = len(dumpInput)
    i = 0
    while i < len(dumpInput):

      if self.output:    
        print "calculating chunks",i+1,"-",min(i+self.cpunum,numberOfRuns),"of",numberOfRuns
      
      if self.debugMode: print "input of this \"map\": ",dumpInput[i:i+self.cpunum]
      
      # create job_server
      job_server = multiprocessing.Pool(processes = self.cpunum)
      
      # map lppWorker on all inputs via job_server (parallelly)
      job_server.map_async(lppWorker,dumpInput[i:i+self.cpunum]).get(9999999)
      
      # close jobserver
      job_server.close()
      job_server.join()
      i += self.cpunum
    
    endtime = time.time()
    if self.output:    
      print "wrote", listlen,"granular snapshots in VTK format"
      print "time needed:",endtime-starttime,"sec"

def lppWorker(input):
  flist = input["filelist"]
  debugMode = input["debugMode"]
  outfileName = input["output"]
  overwrite = input["overwrite"]
  
  flistlen = len(flist)
  # generate name of manyGran
  splitfname = flist[0].rsplit(".")
  if outfileName == "":
    granName = splitfname[len(splitfname)-1]
  elif outfileName.endswith("/"):
    granName = outfileName + splitfname[len(splitfname)-1]
  else:
    granName = outfileName
  
  # if no-overwrite: read timestamp in first line of file
  # if corresponding dump-file does not already exists: add it to 'shortFlist'
  # shortFlist ... list of files to finally be processed by dump, and vtk.
  # elements of flist that are not in shortFlist already exist and will not be
  # converted anew and replaced
  shortFlist = []
  if overwrite == True:
    shortFlist = flist
  else:
    for f in flist:
      try:
        # read time
        ff = open(f)
        ff.readline()
        time = int(ff.readline())
        sys.stdout.flush()
        ff.close()
      except:
        continue
      
      # generate filename from time like in vtk,
      # check if file exists; if yes: do not add to list
      filename,file_bb,file_walls = vtk.generateFilename(granName,[time],0)
      if not os.path.isfile(filename):
        shortFlist.append(f)
  
  # call dump, vtk, manyGran on shortFlist
  try:
    d = dump({"filelist":shortFlist, "debugMode":debugMode})
    v = vtk.vtk(d)
    if debugMode: print "\nfileNums: ",d.fileNums,"\n"
    v.manyGran(granName,fileNos=d.fileNums,output=debugMode)
  except KeyboardInterrupt:
    raise
  
  return 0

def printHelp():
  print "usage: pizza [options] dump.example\n where dump.example is a filename",\
    "or regular expression passing all relevant dump files to pizza."
  print "Important command line options:"
  print "-o fname    : define output file name (default is liggghts + timestep number)"
  print "--chunksize : sets the chunksize, default: 8"
  print "--cpunum    : sets the number of processes to start, default (and maximum)",\
    "is the amout of cpu cores avaliable at your system"
  print "--help      : writes this help message and exits"
  print "--no-overwrite: disables overwriting of already post-processed files."
  print "For details, read README_GRANULAR.txt"

if __name__ == "__main__":
  if len(sys.argv) > 1:
    # parse options
    optlist, args = getopt.gnu_getopt(sys.argv[1:],'o:',['chunksize=','cpunum=','debug','help','quiet','no-overwrite'])
    optdict = dict(optlist)
    if "--help" in optdict:
      printHelp()
    else:
      try:
        lpp(args,**optdict)
      except KeyboardInterrupt:
        print "aborted by user"
      except BaseException, e:
        print "aborting due to errors:", e

    #===========================================================================
    # except:
    #  if sys.exc_type == exceptions.SystemExit:
    #    pass
    #  else: 
    #    print sys.exc_info()
    #===========================================================================
  else:
    printHelp()

  
