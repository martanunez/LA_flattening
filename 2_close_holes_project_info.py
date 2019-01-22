# Close holes corresponding to PVs and LAA. Mark filled holes with a scalar array. Additionally, transfer all scalar arrays from input mesh to output (closed) mesh
# Hole filling done with implementation from https://github.com/cbutakoff/tools/tree/master/FillSurfaceHoles  Related publication: P. Liepa "Filling Holes in Meshes", 2003.
# Hole filling can also be done manually with reMESH (http://remesh.sourceforge.net/)

# Input: mesh with PVs and LAA clipped at their ostia + same mesh without MV
# Output: mesh with holes corresponding to PVs and LAA filled and marked with scalar array. No MV.
# Usage: python 2_close_holes_project_info.py data/mesh_crinkle_clipped.vtk data/mesh_clipped_mitral.vtk data/mesh_clipped_c.vtk

from aux_functions import *
import sys
import os

filename_open = sys.argv[1]
filename_no_mitral = sys.argv[2]
filename_closed = sys.argv[3]

fileroot = os.path.dirname(filename_open)
filename = os.path.basename(filename_open)
filenameroot = os.path.splitext(filename)[0]

if not os.path.exists(filename_open):
    print('ERROR: Input file 1 (LA after PV, LAA clipping) not found')
if not os.path.exists(filename_no_mitral):
    print('ERROR: Input file 2 (LA after PV, LAA, and MV clipping) not found')
if os.path.exists(filename_closed):
    print('WARNING: Closed mesh already exists. Delete it and run again if you want to update it.')
else: # Fill holes
    os.system('./FillSurfaceHoles -i ' + filename_open + ' -o ' + filename_closed)

m_open = readvtk(filename_open)
m_no_mitral = readvtk(filename_no_mitral)
m_closed = readvtk(filename_closed)

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
writevtk(m_final, filename_closed, 'binary')

print('PV and LAA holes have been closed and marked with scalar array <hole> = 1')