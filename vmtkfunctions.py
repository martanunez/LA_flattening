from vmtk import vmtkscripts

def vmtksurfacereader(filename):
    reader = vmtkscripts.vmtkSurfaceReader()
    reader.InputFileName = filename
    reader.Execute()
    return reader.Surface

def vmtkcenterlines(surface, sourcepoints, targetpoints, endpoints=0):
    computer = vmtkscripts.vmtkCenterlines()
    computer.Surface = surface
    computer.SeedSelectorName = 'pointlist'
    computer.SourcePoints = sourcepoints
    computer.TargetPoints = targetpoints
    computer.AppendEndPoints = endpoints
    computer.Execute()
    return computer.Centerlines

def vmtksurfacewriter(polydata, filename):
    # print "Writing",filename
    writer = vmtkscripts.vmtkSurfaceWriter()
    writer.Surface = polydata
    writer.OutputFileName = filename
    writer.Execute()

def vmtkcenterlineresampling(centerline, length=.1):
    resampler = vmtkscripts.vmtkCenterlineResampling()
    resampler.Centerlines = centerline
    resampler.Length = length
    resampler.Execute()
    return resampler.Centerlines

def vmtkcenterlinesmoothing(centerline, iterations=100, factor=0.1):
    smoother = vmtkscripts.vmtkCenterlineSmoothing()
    smoother.Centerlines = centerline
    smoother.NumberOfSmoothingIterations = iterations
    smoother.SmoothingFactor = factor
    smoother.Execute()
    return smoother.Centerlines

def vmtkbranchextractor(centerline):
    extractor = vmtkscripts.vmtkBranchExtractor()
    extractor.Centerlines = centerline
    extractor.RadiusArrayName = 'MaximumInscribedSphereRadius'
    extractor.Execute()
    return extractor.Centerlines

def vmtkcenterlinemerge(centerline, length=.1):
    merger = vmtkscripts.vmtkCenterlineMerge()
    merger.Centerlines = centerline
    merger.Length = length
    merger.RadiusArrayName = 'MaximumInscribedSphereRadius'
    merger.GroupIdsArrayName = 'GroupIds'
    merger.CenterlineIdsArrayName = 'CenterlineIds'
    merger.BlankingArrayName = 'Blanking'
    merger.TractIdsArrayName = 'TractIds'
    merger.Execute()
    return merger.Centerlines

def vmtkcenterlineattributes(centerline):
    computer = vmtkscripts.vmtkCenterlineAttributes()
    computer.Centerlines = centerline
    computer.Execute()
    return computer.Centerlines

def vmtkcenterlinesections(surface, centerline):
    sectioner = vmtkscripts.vmtkCenterlineSections()
    sectioner.Surface = surface
    sectioner.Centerlines = centerline
    sectioner.Execute()
    return sectioner.CenterlineSections