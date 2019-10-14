#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2019 University of Vienna
#
#    This file is part of SHARC.
#
#    SHARC is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SHARC is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    inside the SHARC manual.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************


import sys

def readFile(fileName,option="rb+"):
    """Read File and return text in file
       
       fileName = str, Name of the file to be read (or path+fileName)
       return   = str, Information in fileName
    """

    try:
        data = open(fileName,option)
    except:
        print(" Error reading file: " + fileName)
        sys.exit()

    text = data.read()
    data.close()  
    return text

def writeOutput(fileName,content):
    """ write content to file [fileName]

        fileName = str, Name of the file to be read 
        content  = str, content written to the file
        return   = None
    """
    try:
      OUT = open(fileName,"w")
    except:
      print("Error writing to file: "+fileName)
      sys.exit()
    
    OUT.write(content)
    OUT.close()
    return;
