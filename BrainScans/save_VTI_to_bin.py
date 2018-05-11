import vtk
import sys

reader = vtk.vtkXMLImageDataReader ()
reader.SetFileName('C:\Users\Romin\Documents\materials_100by119by100.vti')
reader.Update()
img = reader.GetOutput()

print "--------------- IMG --------------------------"
print img
print "--------------- img.GetSpacing() -------------"
print img.GetSpacing()
print "--------------- img.GetDimensions() ----------"
print img.GetDimensions()
dim = img.GetDimensions()

print img.GetNumberOfCells()

celldata = img.GetPointData()
print "--------------- celldata.GetNumberOfArrays() -"
print celldata.GetNumberOfArrays()
scalars = celldata.GetArray(0)

print "-----------------SCALARS -----------------"
print scalars

f = open('C:\Users\Romin\Documents\materials_100by119by100.bin','wb+')

status = 'in progress'

import struct

vals = []

counter = 0;

for k in xrange(dim[2]):
    for j in xrange(dim[1]):
        for i in xrange(dim[0]):
            
            counter += 1;
            
            
            if counter%1000 == 0:
                bar_len = 60
                filled_len = int(round(bar_len * counter / float(img.GetNumberOfPoints())))

                percents = round(100.0 * counter / float(img.GetNumberOfPoints()), 1)
                bar = '=' * filled_len + '-' * (bar_len - filled_len)

                sys.stdout.write('[%s] %s%s ...%s %d over %d\r' % (bar, percents, '%', status, counter,img.GetNumberOfPoints()))
                sys.stdout.flush()
                s = struct.pack('f'*len(vals), *vals)
                f.write(s)
                vals = []
            c = scalars.GetTuple(i + dim[0] * ( j + k * dim[1] ) );
            #print c
            vals.append(c[0])


s = struct.pack('f'*len(vals), *vals)
f.write(s)
f.close()

# execfile('C:\\Users\\Romin\\Documents\\save_VTI_to_bin.py')


