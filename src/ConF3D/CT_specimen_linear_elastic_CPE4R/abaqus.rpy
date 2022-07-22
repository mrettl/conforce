# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2020 replay file
# Internal Version: 2019_09_13-19.49.31 163176
# Run by p1746059 on Thu Jul 14 10:50:36 2022
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=301.087493896484, 
    height=116.910003662109)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
openMdb(
    pathName='D:/Tauscher_Markus/Projekte/Configuartional_forces/Version_003/__doku_CF3D/cf3d/src/ConF3D/CT_specimen_linear_elastic/CT_specimen_linear_elastic.cae')
#: The model database "D:\Tauscher_Markus\Projekte\Configuartional_forces\Version_003\__doku_CF3D\cf3d\src\ConF3D\CT_specimen_linear_elastic\CT_specimen_linear_elastic.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
p = mdb.models['CT_specimen'].parts['CT_specimen']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
