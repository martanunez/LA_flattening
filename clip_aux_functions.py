# Code adapted from https://github.com/catactg/SUM
# Related publication describing semi-automatic PV clipping algorithm: "Benchmark for Algorithms Segmenting the Left Atrium From 3D CT and MRI Datasets",
# Catalina Tobon-Gomez et al., IEEE transactions on medical imaging, 2015.

from aux_functions import *
from vmtkfunctions import *
import seedselector

def seed_interactor(surface):
    """Interactor for seed selection. Needs VMTK"""
    computer = seedselector.vmtkPickPointSeedSelector()
    computer.SetSurface(surface)
    computer.Execute()
    return computer.GetSourceSeedIds()

def select_seeds(surface, labels, surfacefileout, vis=0, laa=0):
    """Select 4 seeds (1 per vein) and a 5th one if laa = 1"""
    if laa ==1:
        labelsrange = [76.0, 77.0, 79.0, 78.0, 36.0]
        nseeds = 5
    else:
        # for each PV
        labelsrange = [76.0, 77.0, 79.0, 78.0]
        nseeds = 4

    seeds = seed_interactor(surface)
    # create the pointset
    newpoints = vtk.vtkPoints()
    newvertices = vtk.vtkCellArray()

    # create array on seeds with ground truth (GT) labels
    gtlabels_array = vtk.vtkDoubleArray()
    gtlabels_array.SetName(labels)

    if not seeds.GetNumberOfIds() == nseeds:
        print('You should select extactly', nseeds, ' seeds. Try again!')
        seeds = seed_interactor(surface)

    for s in range(seeds.GetNumberOfIds()):
        branchlabel = labelsrange[s]
        point = surface.GetPoint(seeds.GetId(s))
        pid = newpoints.InsertNextPoint(point)
        gtlabels_array.InsertNextValue(branchlabel)
        # Create the topology of the point (a vertex)
        newvertices.InsertNextCell(1)
        newvertices.InsertCellPoint(pid)
    pointspd = vtk.vtkPolyData()
    pointspd.SetPoints(newpoints)
    pointspd.SetVerts(newvertices)
    pointspd.GetPointData().AddArray(gtlabels_array)

    if vis==1:
        pointsgplyh = generateglyph(pointspd)
        visualise_default(pointsgplyh, surface, 'seeds', labels, 36, 79)
    writevtp(pointspd, surfacefileout)

def visualise_default(surface, ref, case, arrayname, mini, maxi):
    """Visualise surface with a default parameters"""
    #Create a lookup table to map cell data to colors
    # print "Colormap from ", mini, "to", maxi
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(255)
    lut.SetValueRange(0, 255)

    # qualitative data from colorbrewer  --> matching qualitative colormap of Paraview
    lut.SetTableValue(0, 0, 0, 0, 1)  #Black
    lut.SetTableValue(mini, 1, 1, 1, 1)   #white
    lut.SetTableValue(mini+1, 77/255.,175/255., 74/255., 1)   # green
    lut.SetTableValue(maxi-3, 152/255.,78/255.,163/255., 1)  # purple
    lut.SetTableValue(maxi-2, 255/255.,127/255., 0., 1)  # orange
    lut.SetTableValue(maxi-1, 55/255., 126/255., 184/255., 1)  # blue
    lut.SetTableValue(maxi, 166/255., 86/255., 40/255., 1)  # brown
    lut.Build()

    # create a text actor
    txt = vtk.vtkTextActor()
    txt.SetInput(case)
    txtprop=txt.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(18)
    txtprop.SetColor(0, 0, 0)
    txt.SetDisplayPosition(20, 30)

    # create a rendering window, renderer, and renderwindowinteractor
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    # for GIMIAS interaction style
    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)
    iren.SetRenderWindow(renWin)

    # surface mapper and actor
    surfacemapper = vtk.vtkPolyDataMapper()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        surfacemapper.SetInputData(surface)
    else:
        surfacemapper.SetInput(surface)
    surfacemapper.SetScalarModeToUsePointFieldData()
    surfacemapper.SelectColorArray(arrayname)
    surfacemapper.SetLookupTable(lut)
    surfacemapper.SetScalarRange(0,255)
    surfaceactor = vtk.vtkActor()
    # surfaceactor.GetProperty().SetOpacity(0)
    # surfaceactor.GetProperty().SetColor(1, 1, 1)
    surfaceactor.SetMapper(surfacemapper)

    # refsurface mapper and actor
    refmapper = vtk.vtkPolyDataMapper()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        refmapper.SetInputData(ref)
    else:
        refmapper.SetInput(ref)
    refmapper.SetScalarModeToUsePointFieldData()
    refmapper.SelectColorArray(arrayname)
    refmapper.SetLookupTable(lut)
    refmapper.SetScalarRange(0,255)
    refactor = vtk.vtkActor()
    refactor.GetProperty().SetOpacity(0.5)
    # refactor.GetProperty().SetColor(1, 1, 1)
    refactor.SetMapper(refmapper)

    # assign actors to the renderer
    ren.AddActor(refactor)
    ren.AddActor(surfaceactor)
    ren.AddActor(txt)

    # set the background and size; zoom in; and render
    ren.SetBackground(1, 1, 1)
    renWin.SetSize(1280, 960)
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1)

    # before
    # print("before", ren.GetActiveCamera().GetViewUp())

    # enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()

    outcam = ren.GetActiveCamera()
    # print("after", outcam.GetViewUp())

def visualise_color(surface, ref, case):
    """Visualise surface in solid color and 'ref' in trasparent"""
    # create a text actor
    txt = vtk.vtkTextActor()
    txt.SetInput(case)
    txtprop=txt.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(18)
    txtprop.SetColor(0, 0, 0)
    txt.SetDisplayPosition(20, 30)

    # create a rendering window, renderer, and renderwindowinteractor
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    # for GIMIAS interaction style
    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)
    iren.SetRenderWindow(renWin)

    # surface mapper and actor
    surfacemapper = vtk.vtkPolyDataMapper()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        surfacemapper.SetInputData(surface)
    else:
        surfacemapper.SetInput(surface)
    surfacemapper.SetScalarModeToUsePointFieldData()
    surfaceactor = vtk.vtkActor()
    # surfaceactor.GetProperty().SetOpacity(0)
    surfaceactor.GetProperty().SetColor(288/255, 26/255, 28/255)
    surfaceactor.SetMapper(surfacemapper)

    # refsurface mapper and actor
    refmapper = vtk.vtkPolyDataMapper()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        refmapper.SetInputData(ref)
    else:
        refmapper.SetInput(ref)
    refmapper.SetScalarModeToUsePointFieldData()

    refactor = vtk.vtkActor()
    refactor.GetProperty().SetOpacity(0.5)
    refactor.GetProperty().SetColor(1, 1, 1)
    refactor.SetMapper(refmapper)

    # assign actors to the renderer
    # ren.AddActor(refactor)
    ren.AddActor(surfaceactor)
    ren.AddActor(refactor)
    ren.AddActor(txt)

    # set the background and size; zoom in; and render
    ren.SetBackground(1, 1, 1)
    renWin.SetSize(800, 800)
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1)

    # enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()

def getregionslabels():
    """Return dictionary linking regionids to anatomical locations."""
    regionslabels = {'body': 36,
                     'laa': 37,
                     'pv2': 76,
                     'pv1': 77,
                     'pv3': 78,
                     'pv4': 79}
    return regionslabels

def create_autolabels(surface, ref, arrayname, value):
    """Create autolabels scalar array (mark PVs using branch labels) and add it to surface """
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(surface)
    locator.BuildLocator()

    array = surface.GetPointData().GetArray(arrayname)
    for i in range(ref.GetNumberOfPoints()):
        point = ref.GetPoint(i)
        closestpoint_id = locator.FindClosestPoint(point)
        array.SetValue(closestpoint_id, value)
    return surface

def centroidofcentroids(edges):
    # compute centroids of each edge
    # find average point
    acumvector = [0,0,0]
    rn = countregions(edges)
    # print "found",rn,'edges'
    for r in range(rn):
        oneedge = extractconnectedregion(edges,r)
        onecentroid = pointset_centreofmass(oneedge)
        acumvector = acumvectors(acumvector,onecentroid)
        # print acumvector
    finalcentroid = dividevector(acumvector,rn)
    return finalcentroid

def pv_LAA_centerlines(inputfile, seedsfile, outfile, pvends=1):
    """ Create 5 pairs of centerlines, each one starting from each PV (or LAA) seed and going to the 2 opposite
    (other side) PVs"""

    # surface = vmtksurfacereader(inputfile)
    surface = readvtk(inputfile)
    points = np.loadtxt(seedsfile, delimiter=',').tolist()

    print('Processing RSPV seed:')
    cl1 = vmtkcenterlines(surface, points[0], points[2] + points[3], pvends)
    print('\n \nProcessing RIPV seed:')
    cl2 = vmtkcenterlines(surface, points[1], points[2] + points[3], pvends)
    print('\n \nProcessing LIPV seed:')
    cl3 = vmtkcenterlines(surface, points[2], points[0] + points[1], pvends)
    print('\n \nProcessing LSPV seed:')
    cl4 = vmtkcenterlines(surface, points[3], points[0] + points[1], pvends)
    print('\n \nProcessing LAA seed:')
    cl5 = vmtkcenterlines(surface, points[4], points[0] + points[1], pvends)

    writevtp(cl1, outfile + 'clraw21.vtp')
    writevtp(cl2, outfile + 'clraw22.vtp')
    writevtp(cl3, outfile + 'clraw23.vtp')
    writevtp(cl4, outfile + 'clraw24.vtp')
    writevtp(cl5, outfile + 'clraw25.vtp')

def intersectwithline(surface, p1, p2):
    """Given surface and line defined by 2 points (p1,p2), return insersecting points"""
    tree = vtk.vtkOBBTree()
    tree.SetDataSet(surface)
    tree.BuildLocator()

    intersectPoints = vtk.vtkPoints()
    intersectCells = vtk.vtkIdList()

    tolerance=1.e-3
    tree.SetTolerance(tolerance)
    tree.IntersectWithLine(p1, p2, intersectPoints, intersectCells)
    return intersectPoints

def furthest_point_to_polydata(pointset,refpoint):
    """Given set of points and ref point, select furthest point using euclidean dist"""
    refdist = 0
    for i in range(pointset.GetNumberOfPoints()):
        dist = euclideandistance(pointset.GetPoint(i),refpoint)
        if dist > refdist:
            refdist = dist
            selectedpointid = i
    return pointset.GetPoint(selectedpointid)

def computelengthalongvector(polydata, refpoint, vector):
    # polydata should be a closed surface

    # intersect with line
    point1 = refpoint
    point2 = sumvectors(refpoint,1000,vector) # far away point
    intersectpoints = intersectwithline(polydata,point1,point2)
    furthestpoint1 = furthest_point_to_polydata(intersectpoints,refpoint)

    # intersect with line the other way
    point1 = refpoint
    point2 = sumvectors(refpoint,-1000,vector) # far away point
    intersectpoints = intersectwithline(polydata,point1,point2)
    furthestpoint2 = furthest_point_to_polydata(intersectpoints,furthestpoint1)
    length = euclideandistance(furthestpoint1,furthestpoint2)
    return length

def clip_vein(surface,cl,clippointid):
    """Clip the vein at clippoint"""
    clippoint0 = cl.GetPoint(clippointid)
    clipnormal = (np.array(cl.GetPoint(clippointid+1)) - np.array(cl.GetPoint(clippointid)))
    possvein = planeclip(surface, clippoint0, clipnormal)
    vein = extractclosestpointregion(possvein,clippoint0)
    return vein

def skippoints(polydata, nskippoints):
    """Generate a single cell line from points in idlist."""
    # derive number of nodes
    numberofnodes = polydata.GetNumberOfPoints() - nskippoints

    # define points and line
    points = vtk.vtkPoints()
    polyline = vtk.vtkPolyLine()
    polyline.GetPointIds().SetNumberOfIds(numberofnodes)

    # assign id and x,y,z coordinates
    for i in range(nskippoints,polydata.GetNumberOfPoints()):
        pointid = i - nskippoints
        polyline.GetPointIds().SetId(pointid,pointid)
        point = polydata.GetPoint(i)
        points.InsertNextPoint(point)

    # define cell
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyline)

    # add to polydata
    polyout = vtk.vtkPolyData()
    polyout.SetPoints(points)
    polyout.SetLines(cells)
    if not vtk.vtkVersion.GetVTKMajorVersion() > 5:
        polyout.Update()
    return polyout

def clip_veins_sections_and_LAA(inputfile, sufixfile, clspacing, maxslope, skippointsfactor, highslope, bumpcriterion):
    """ We wish to clip the vein as close to the body as possible without
    including parts of the body or other veins. 'Trial' clips are
    obtained using vmtkcenterlinesections, which creates for each point
    on the centerline a section perpendicular to the centerline and
    provides measures of the section such as maximum diameter.
    When the series of sections enter the atrium body, the maximum
    diameter increases significantly. To quantify the change in max
    diameter between one section and the next in terms of centerline
    spacing, we define 'slope'. When this slope exceeds a certain
    threshold, we assume to have entered the body. The clippoint is
    defined as the centerline point corresponding to the last section of
    the vein before entering the body. """

    surface = vmtksurfacereader(inputfile)  # atrium surface mesh
    # surface = readvtk(inputfile)  # atrium surface mesh

    # creating array to hold new autolabels
    branchlabel= [0, 77, 76, 78, 79, 37]    # it must be 37 for the LAA

    branch_array = vtk.vtkDoubleArray()
    branch_array.SetName('autolabels')
    branch_array.SetNumberOfTuples(surface.GetNumberOfPoints())
    surface.GetPointData().AddArray(branch_array)

    # initialize with bodylabel
    for i in range(surface.GetNumberOfPoints()):
        branch_array.SetValue(i, round(36))

    #for k in range(1, 5):
    for k in range(1, 6):
        print("branchlabel", branchlabel[k])
        cl = readvtp(sufixfile + 'clraw2' + str(k) + '.vtp')  # 2 means with endpoints
        cl = vmtkcenterlineresampling(cl, clspacing)
        cl = vmtkcenterlinesmoothing(cl)
        cl = vmtkbranchextractor(cl)
        writevtp(cl, sufixfile + 'clbranch' + str(k) + '.vtp')

        #-----------------------------------------------------------------------
        # BUG FIX
        # for some anatomies, the branch extractor gives different
        # number of branches for both centerlines.
        # Calling cellthreshold seems to fix this, instead of previous
        # fix:
        # groupids = [0, 1, 2, 0, 1, 3]
        # for i in range(cl.GetNumberOfCells()):
        #     cl.GetCellData().GetArray('GroupIds').SetValue(i, groupids[i])
        #-----------------------------------------------------------------------
        cl = cellthreshold(cl, 'GroupIds', 0, 0)
        cl = vmtkcenterlinemerge(cl)

        # cl = extractcells(cl, [0])  # vein's centerline    BUG in Windows. If not commented, cl is empty and python crashes without error message

        # original cl for clipping
        cl = vmtkcenterlineresampling(cl, clspacing)
        cl = vmtkcenterlineattributes(cl)
        writevtp(cl, sufixfile +'clvein' + str(k) + '.vtp')

        # removing endpoints that cause failure in vmtkcenterlinesections
        # cl only used in sections
        # print cl.GetNumberOfPoints(),skippointsfactor
        nskippoints = round(skippointsfactor*cl.GetNumberOfPoints())
        # print "skipping ", nskippoints, " points from ", cl.GetNumberOfPoints()
        cl = skippoints(cl, int(nskippoints))
        cl = vmtkcenterlineresampling(cl, clspacing)
        cl = vmtkcenterlineattributes(cl)
        #-----------------------------------------------------------------------
        # Tiral clips
        # The algorithm is constructed such that the clippoint can not be the
        # first or last point of the centerline. This property is required for
        # calculating clipnormal as defined below.
        #-----------------------------------------------------------------------
        # print "pre sections"
        sections = vmtkcenterlinesections(surface, cl)
        # print "post sections"
        vmtksurfacewriter(sections, sufixfile + 'clsection' + str(k) + '.vtp')

        closedarray = sections.GetCellData().GetArray('CenterlineSectionClosed')
        maxsizearray = (sections.GetCellData().GetArray('CenterlineSectionMaxSize'))

        highcount = 0
        nbumpcriterion  = round(bumpcriterion*cl.GetNumberOfPoints())

        print("bumps ", nbumpcriterion, " points from ", cl.GetNumberOfPoints())
        for i in range(1, sections.GetNumberOfCells()):

            # Skip sections that are preceded by open sections and skip first
            # 5 centerline points. This to avoid complications near the end of a
            # vein (holes in the surface mesh or sudden changes in centerline
            # direction might otherwise lead to exceeding the threshold far
            # from the atrium body).
            if closedarray.GetValue(i-1):
                # changed to "signed" difference to account for veins w
                # multiple outlets. This created a thin to wide vein change

                slope = (maxsizearray.GetValue(i) - maxsizearray.GetValue(i-1))/ clspacing

                if slope > highslope:
                    highcount += 1
                else:
                    highcount = 0

                # print i, slope, highcount

                if slope > maxslope:
                    break
                elif slope > highslope and highcount == nbumpcriterion:
                    break
                else:
                    pass

        if highcount == 0:
            clippointid = i - 1
        else:
            clippointid = i - (highcount)

        # print clippointid
        np.savetxt(sufixfile + 'clippointid' + str(k) + '.csv', np.array([clippointid])+int(nskippoints), fmt='%i')

        #-----------------------------------------------------------------------
        # Prepare output
        #
        # Tried clipping to make nice ostia clips, but it's not so trivial.
        # Back to good old transfer labels
        #-----------------------------------------------------------------------
        vein = clip_vein(surface, cl, clippointid)
        surface = create_autolabels(surface, vein, 'autolabels', round(branchlabel[k]))
        # surface = transfer_array(surface, vein, 'autolabels', 'autolabels')
        # visualise_color(vein,surface,'vein' + str(k))
    writevtp(surface, sufixfile + 'autolabels.vtp')

def clip_vein_endpoint_and_LAA_save_planes(surface, ifile_sufix, targetdistance, specialvein=0, specialdist=0):
    """Clip vein the targetdistance away from the body. Clip also the LAA at specialdist.
    Return the clip planes, for each plane: point + normal
    in a numpy matrix. First row = 1st point (x,y,z), Second row = 1st normal (x,y,z). Then continue with the rest of PVs and LAA
    """
    clip_planes = np.zeros((10,3))
    regionslabels = getregionslabels()

    # extract the body from the surface
    body = pointthreshold(surface, 'autolabels', regionslabels['body'], regionslabels['body'], 1)
    body = extractlargestregion(body)

    # initialize appender with the body
    appender = vtk.vtkAppendPolyData()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        appender.AddInputData(body)
    else:
        appender.AddInput(body)
    originaldist = targetdistance
    for k in range(1, 6):
        if k == 5:
            index = 'laa'
        else:
            index = 'pv' + str(k)
        # extract vein
        # excluding some points (alloff=0)
        # to avoid overlapping edges after appending
        vein = pointthreshold(surface, 'autolabels', regionslabels[index], regionslabels[index], 0)

        # load the centreline and the clipoint
        cl = readvtp(ifile_sufix + 'clvein' + str(k) + '.vtp')
        clippointid = int(np.loadtxt(ifile_sufix + 'clippointid' + str(k) + '.csv'))

        clippoint0 = cl.GetPoint(clippointid)
        clipnormal = (np.array(cl.GetPoint(clippointid + 1)) - np.array(cl.GetPoint(clippointid )))

        abscissasarray = cl.GetPointData().GetArray('Abscissas')
        startabscissa = abscissasarray.GetValue(clippointid)
        currentabscissa = 0
        currentid = clippointid

        # if different distance for 1 vein
        if specialvein > 0:
            if regionslabels[index] == specialvein:
                targetdistance = specialdist
            else:
                targetdistance = originaldist

        # find clip point
        while ((currentabscissa < targetdistance) and (currentabscissa >= 0) and (currentid >= 0)):
            currentid -= 1
            currentabscissa = startabscissa - abscissasarray.GetValue(currentid)

        if currentid > 0:
            currentid = currentid + 1
        else:
            # vein ended before target distance
            # then clip 2 mm before end of centreline from end point
            currentid = 4

        # clip and append
        clippoint1 = cl.GetPoint(currentid)
        clippedvein = planeclip(vein, clippoint1, clipnormal, 0)

        clip_planes[2*(k-1), 0:3] = clippoint1
        clip_planes[2*(k-1) +1, 0:3] = clipnormal

        # keep region closest to ostium point
        clippedvein = extractclosestpointregion(clippedvein, clippoint0)

        # clip generates new points to make a flat cut. The values may be interpolated.
        # we want all values to rounded to a certain label value.
        clippedvein = roundpointarray(clippedvein, 'autolabels')
        if vtk.vtkVersion.GetVTKMajorVersion() > 5:
            appender.AddInputData(clippedvein)
        else:
            appender.AddInput(clippedvein)

    # collect body + veins
    appender.Update()
    clippedsurface = appender.GetOutput()
    clippedsurface = cleanpolydata(clippedsurface)
    return clippedsurface, clip_planes

def cylinderclip(dataset, point0, point1,normal,radius):
    """Define cylinder. The cylinder is infinite in extent. We therefore have
    to truncate the cylinder using vtkImplicitBoolean in combination with
    2 clipping planes located at point0 and point1. The radius of the
    cylinder is set to be slightly larger than 'maxradius'."""

    rotationaxis = cross([0, 1, 0], normal)
    rotationangle = (180 / math.pi) * angle([0, 1, 0], normal)

    transform = vtk.vtkTransform()
    transform.Translate(point0)
    transform.RotateWXYZ(rotationangle, rotationaxis)
    transform.Inverse()

    cylinder = vtk.vtkCylinder()
    cylinder.SetRadius(radius)
    cylinder.SetTransform(transform)

    plane0 = vtk.vtkPlane()
    plane0.SetOrigin(point0)
    plane0.SetNormal([-x for x in normal])
    plane1 = vtk.vtkPlane()
    plane1.SetOrigin(point1)
    plane1.SetNormal(normal)

    clipfunction = vtk.vtkImplicitBoolean()
    clipfunction.SetOperationTypeToIntersection()
    clipfunction.AddFunction(cylinder)
    clipfunction.AddFunction(plane0)
    clipfunction.AddFunction(plane1)

    clipper = vtk.vtkClipPolyData()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        clipper.SetInputData(dataset)
    else:
        clipper.SetInput(dataset)
    clipper.SetClipFunction(clipfunction)
    clipper.Update()
    return extractlargestregion(clipper.GetOutput())

def find_mitral_cylinder_pvs(surface, arrayname, outfile, scale=0.4, w=[0.7,0.15,0.15], vis=0):
    """Compute local coordinate system based on the body centroid and PVs centroid.
    The 3 axes are weighted as in w. The resulting vector is used to clip surface
    scale * radius away from the body centroid."""

    # extract body
    startedges = extractboundaryedge(surface)
    if startedges.GetNumberOfPoints() > 0:
        surfacefilled = fillholes(surface, 1000)
    else:
        surfacefilled = surface

    body = pointthreshold(surface, arrayname, 36.0, 36.0)
    bodycom = pointset_centreofmass(body)

    # average of left ostia to average of right ostia
    ostia = pointthreshold(surfacefilled, arrayname, 78.0, 79.0)
    edges = extractboundaryedge(ostia)

    # print "pvs"
    leftcentroid = centroidofcentroids(edges)
    # print 'leftcentroid',leftcentroid

    ostia = pointthreshold(surfacefilled, arrayname, 76.0, 77.0)
    edges = extractboundaryedge(ostia)

    rightcentroid = centroidofcentroids(edges)

    # final pvscom average of left and right
    pvscom = acumvectors(leftcentroid, rightcentroid)
    pvscom = dividevector(pvscom, 2)

    # NOW AXES
    # Axis 1: Pvs com to body com
    pvdir = subtractvectors(bodycom, pvscom)
    pvdirn = normalizevector(pvdir)

    # Axis 2: normal to Pvs axis
    ostiadir1 = subtractvectors(leftcentroid, rightcentroid)
    ostiadirn = normalizevector(ostiadir1)

    ostiacross = cross(pvdirn, ostiadirn)
    ostiacrossn = normalizevector(ostiacross)

    # Axis 3: normal to axis 1 and 2
    pvcross = cross(ostiacrossn, pvdirn)
    pvcrossn = normalizevector(pvcross)

    # thought of using for weighting but defualt values seem all right
    bodylength= computelengthalongvector(body, bodycom, pvdirn)
    measurepoint = sumvectors(bodycom, scale*bodylength, pvdirn)
    bodythick = computelengthalongvector(body, measurepoint, ostiacrossn)
    bodywidth = computelengthalongvector(body, measurepoint, pvcrossn)

    # print('length', bodylength, 'widht', bodywidth, 'thickness', bodythick)

    pvdirnw = multiplyvector(pvdirn, w[0])
    ostiadirnw = multiplyvector(pvcrossn, w[1])
    ostiacrossnw = multiplyvector(ostiacrossn, w[2])

    plusvector = acumvectors(pvdirnw, ostiacrossnw)
    plusvector = acumvectors(plusvector, ostiadirnw)
    plusvectorn = normalizevector(plusvector)

    # clippoint with length vector
    # in very small bodies, modify scale
    if bodylength/bodythick < 1.5:
        print("short body")
        scale = 0.45
    clippoint = sumvectors(bodycom, scale*bodylength, plusvectorn)

    # current cut
    slicepv = cutdataset(surface, clippoint, plusvectorn)
    # if > 1 region --> clipping LAA as well.
    nr = countregions(slicepv)
    if nr > 1:
        print("recalculating clip point")
        # keep cut closest to point
        slicepv = extractclosestpointregion(slicepv, clippoint)
        # recompute center
        clippoint = pointset_centreofmass(slicepv)

    vectordown = multiplyvector(plusvectorn, 10)
    pointdown = acumvectors(vectordown, clippoint)

    finalbody = cylinderclip(surface, clippoint, pointdown, plusvectorn, bodythick)
    # print 'point0 = ', clippoint
    # print 'point1 = ', pointdown
    # print 'normal = ', plusvectorn
    # print 'radius = ', bodythick

    if vis == 1:
        # Visualisation of coordinate system
        pvscompd = point2vertexglyph(pvscom)
        pvscompdg = generateglyph(pvscompd)

        plotpoint = sumvectors(pvscom, w[0]*bodylength, pvdirn)
        bodyaxis = linesource(plotpoint, pvscom)
        plotpoint = sumvectors(pvscom, w[1]*bodylength, pvcrossn)
        ostiaaxis = linesource(plotpoint, pvscom)
        plotpoint = sumvectors(pvscom, w[2]*bodylength, ostiacrossn)
        crossaxis = linesource(plotpoint, pvscom)

        allaxis = append(pvscompdg, ostiaaxis)
        allaxis = append(allaxis, bodyaxis)
        allaxis = append(allaxis, crossaxis)
        allaxis = append(allaxis, slicepv)

        writevtp(allaxis, outfile + '_axes.vtp')
        visualise_color(allaxis, surfacefilled, 'plus')
    # writevtp(finalbody, outfile + '.vtp')
    return finalbody
