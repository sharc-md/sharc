#!/usr/remote/bin/python -u

# Crude script for the conversion of
# orca output files to frequency files
# for molden. Only the last occurring
# instances of coordinates, frequencies
# and normal modes are used.

from sys import argv

# Check arguments
if len( argv ) != 2:
    print "Usage: orca2molden <orca-outfile>"
    print ""
    print "Convert orca output to molden file for"
    print "normal mode visualisation."
    exit()

name, orca_file = argv

try:
    lines = open( orca_file, 'r' ).readlines()
except IOError:
    print "Could not open %s." %orca_file
    exit()

# check if file is sucessfully completed orca file:
is_orca = False
finished = False
if lines[2].strip() == "* O   R   C   A *":
    is_orca = True
if lines[-2].strip() == "****ORCA TERMINATED NORMALLY****":
    finished = True

if not is_orca:
    print "File %s is not in orca output format (probably)!" % orca_file
    exit()
elif is_orca and not finished:
    print "The job either has crashed or has not finished yet."
    exit()
elif is_orca and finished:
    print "Reading data from file %s..." % orca_file

# set falgs
read = False
geom = False
freq = False
intensity=False
nmode = False
empty = 0
counter = 0

# parse input
for line in lines:
    if line.startswith( "CARTESIAN COORDINATES (A.U.)" ):
        read = True
        geom = True
        coords = []
    elif line.startswith( "VIBRATIONAL FREQUENCIES" ):
        read = True
        freq = True
        freqs = []
        empty = 0
    elif line.startswith( "Thus, these vectors are normalized but" ):
        read = True
        nmode = True
        modes = []
    elif line.startswith( "IR SPECTRUM" ):
        read=True
        intensity=True
        intdict={}
    elif read == True and geom == True:
        line = line.strip()
        if line.startswith( "----" ) or line.startswith( "NO" ):
            continue
        if line:
            line = line.split()
            coords.append( "%s %s %s %s" %( line[1], line[5], line[6], line[7] ) )
        else:
            read = False
            geom = False
    elif read == True and freq == True:
        line = line.strip()
        if line.startswith( "----" ):
           continue
        elif not line and empty == 0:
            empty += 1
            continue
        elif line:
            line = line.split()
            freqs.append( line[1] )
        elif not line and empty > 0:
            read = False
            freq = False
            empty = 0
        else:
            read = False
            freq = False
    elif read == True and nmode == True:
        n_atoms = len( coords )
        cart = n_atoms * 3
        line = line.strip()
        if not line:
            continue
        elif line.startswith( "----" ):
            counter = 0
            read = False
            nmode = False
        elif line:
            counter += 1
            if ( counter+cart )  % (cart+1) == 0:
                # remove normal mode headers
                pass
            else:
                line = line.split()[1:]
                modes.append( line )
        else:
            counter = 0
            read = False
            nmode = False
    elif read == True and intensity == True:
        line = line.strip()
        if line=='':
            continue
        if line.startswith( "----" ):
           continue
        if 'Mode' in line:
           continue
        if 'The first' in line:
           intensity=False
           read=False
           continue
        l=line.replace(':','').strip().split()
        intdict[int(l[0])]=float(l[2])


# check coords
if len( coords ) > 0:
    print "Found system of %d atoms." % len( coords )
else:
    print "No cartesian coordinates found."
    exit()
# check freqs
if len( freqs ) > 0:
    print "Found %d vibrational frequencies." %len( freqs )
else:
    print "No vibrational frequencies found."
# check modes
if len( modes ) == 0:
    print "No normal modes found."
    exit()
# check intensities
if len(intdict) == 0:
    print 'No intensities available.'
else:
    print 'IR intensities available.'

# combine modes
n_atoms = len( coords )
cart = n_atoms * 3
for i in xrange( cart, len( modes ) ):
    modes[i%cart] = modes[i%cart] + modes[i]

# sort modes
n_modes = len( modes[0] )
print "Found %d normal modes." % n_modes
modes = modes[:n_modes]
cmodes = [ [] for i in xrange( n_modes ) ]
for line in modes:
    for i in xrange( n_modes ):
        cmodes[i].append( line[i] )
modes = cmodes

# check for consistency:
if n_modes != len( freqs ):
    print "Mismatch between frequencies and vibrational normal modes detected."
    exit()

# generate molden file
out_file = orca_file + '.molden'
out = open( out_file, 'w' )
out.write( "[MOLDEN FORMAT]\n" )
# write frequencies
out.write( "[FREQ]\n" )
for freq in freqs:
    out.write( freq+'\n' )
# write coordinates block (A.U.)
out.write( "[FR-COORD]\n" )
for coord in coords:
    out.write( coord+'\n' )
# write normal modes:
out.write( "[FR-NORM-COORD]\n" )
for i in xrange( n_modes ):
    out.write( "vibration %d\n" %(i+1) )
    for j in xrange( len( modes[i] ) ):
        out.write( modes[i][j]+' ' )
        if (j+1)%3 == 0:
            out.write( '\n' )
out.write('[INT]\n')
for i in xrange(n_modes):
    if i in intdict:
        f=intdict[i]
    else:
        f=0.
    out.write('%16.9f\n' % f)
out.close()
print "Molden output written to %s" % out_file






# kate: indent-width 4
