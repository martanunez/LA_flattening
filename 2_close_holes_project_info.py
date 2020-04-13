# Close holes corresponding to PVs and LAA. Mark filled holes with a scalar array. Additionally, transfer all scalar arrays from input mesh to output (closed) mesh
# Hole filling done with implementation from https://github.com/cbutakoff/tools/tree/master/FillSurfaceHoles  Related publication: P. Liepa "Filling Holes in Meshes", 2003.
# Hole filling can also be done manually with reMESH (http://remesh.sourceforge.net/)

# Input: mesh with PVs and LAA clipped at their ostia + same mesh without MV
# Output: mesh with holes corresponding to PVs and LAA filled and marked with scalar array. No MV.
# Usage: python 2_close_holes_project_info.py --meshfile_open data/mesh_crinkle_clipped.vtk --meshfile_open_no_mitral  data/mesh_clipped_mitral.vtk --meshfile_closed data/mesh_clipped_c.vtk

from aux_functions import *
import sys
import os
import argparse
from sys import platform

parser = argparse.ArgumentParser()
parser.add_argument('--meshfile_open', type=str, metavar='PATH', help='path to input mesh with clipped PVs and LAA')
parser.add_argument('--meshfile_open_no_mitral', type=str, metavar='PATH', help='path to input mesh with additional MV clip')
parser.add_argument('--meshfile_closed', type=str, metavar='PATH', help='path to output mesh, i.e. with filled holes')
args = parser.parse_args()

fileroot = os.path.dirname(args.meshfile_open)
filename = os.path.basename(args.meshfile_open)
filenameroot = os.path.splitext(filename)[0]

if not os.path.exists(args.meshfile_open):
    sys.exit('ERROR: Input file 1 (LA after PV, LAA clipping) not found')
if not os.path.exists(args.meshfile_open_no_mitral):
    sys.exit('ERROR: Input file 2 (LA after PV, LAA, and MV clipping) not found')
if os.path.exists(args.meshfile_closed):
    print('WARNING: Closed mesh already exists. Delete it and run again if you want to update it.')
else:  # Fill holes
    #os.system('./FillSurfaceHoles -i ' + args.meshfile_open + ' -o ' + args.meshfile_closed)
    if platform == "linux" or platform == "linux2":
        os.system('./FillSurfaceHoles -i ' + args.meshfile_open + ' -o ' + args.meshfile_closed)
    elif platform == "win32":
        os.system('FillSurfaceHoles_Windows\FillSurfaceHoles.exe -i ' + args.meshfile_open + ' -o ' + args.meshfile_closed + ' -smooth none')   # default smooth cotangent (and edglen) fails when using this binary
    else:
        sys.exit('Unknown operating system. Holes cannot be filled automatically. Fill holes manually and save file as ', args.meshfile_closed, '. Then run again this script to proyect scalar arrays from initial mesh if necessary.')

m_open = readvtk(args.meshfile_open)
m_no_mitral = readvtk(args.meshfile_open_no_mitral)
m_closed = readvtk(args.meshfile_closed)

print('Projecting information... ')
transfer_all_scalar_arrays(m_open, m_closed)

# Mark filled holes. Common points (close enough, not added during hole filling) will we marked with scalar array
array_labelsA = np.zeros(m_closed.GetNumberOfPoints())
locator = vtk.vtkPointLocator()
locator.SetDataSet(m_open)
locator.BuildLocator()
for p in range(m_closed.GetNumberOfPoints()):
    point = m_closed.GetPoint(p)
    closestpoint_id = locator.FindClosestPoint(point)
    dist = euclideandistance(point, m_open.GetPoint(closestpoint_id))
    if dist > 0.05:   # empirical distance
        array_labelsA[p] = 1
newarray = numpy_to_vtk(array_labelsA)
newarray.SetName('hole')
m_closed.GetPointData().AddArray(newarray)

# Mark MV using m_open and m_no_mitral
array_labelsB = np.zeros(m_open.GetNumberOfPoints())
locator = vtk.vtkPointLocator()
locator.SetDataSet(m_no_mitral)
locator.BuildLocator()
for p in range(m_open.GetNumberOfPoints()):
    point = m_open.GetPoint(p)
    closestpoint_id = locator.FindClosestPoint(point)
    dist = euclideandistance(point, m_no_mitral.GetPoint(closestpoint_id))
    if dist > 0.05:   # empirical distance
        array_labelsB[p] = 1
newarray = numpy_to_vtk(array_labelsB)
newarray.SetName('mv')
m_open.GetPointData().AddArray(newarray)

transfer_array(m_open, m_closed, 'mv', 'mv')
transfer_array(m_open, m_closed, 'autolabels', 'autolabels')
m_final = pointthreshold(m_closed, 'mv', 0, 0)
m_final.GetPointData().RemoveArray('mv')
writevtk(m_final, args.meshfile_closed, 'binary')

print('PV and LAA holes have been closed and marked with scalar array <hole> = 1')
