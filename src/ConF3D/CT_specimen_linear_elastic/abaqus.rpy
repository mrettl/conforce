# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2020 replay file
# Internal Version: 2019_09_13-19.49.31 163176
# Run by p1746059 on Thu Jul 21 17:42:22 2022
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
o1 = session.openOdb(
    name='D:/Tauscher_Markus/Projekte/Configuartional_forces/Version_003/__doku_CF3D/cf3d/src/ConF3D/CT_specimen_linear_elastic_CPE4/CT_specimen_CPE4.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#* OdbError: Cannot open file 
#* D:/Tauscher_Markus/Projekte/Configuartional_forces/Version_003/__doku_CF3D/cf3d/src/ConF3D/CT_specimen_linear_elastic_CPE4/CT_specimen_CPE4.odb. 
#* *** ERROR: No such file: 
#* D:/Tauscher_Markus/Projekte/Configuartional_forces/Version_003/__doku_CF3D/cf3d/src/ConF3D/CT_specimen_linear_elastic_CPE4/CT_specimen_CPE4.odb.
o1 = session.openOdb(
    name='D:/Tauscher_Markus/Projekte/Configuartional_forces/Version_003/__doku_CF3D/cf3d/src/ConF3D/CT_specimen_linear_elastic/CT_specimen_CPE4.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: D:/Tauscher_Markus/Projekte/Configuartional_forces/Version_003/__doku_CF3D/cf3d/src/ConF3D/CT_specimen_linear_elastic/CT_specimen_CPE4.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       5
#: Number of Node Sets:          37
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
odb = session.odbs['D:/Tauscher_Markus/Projekte/Configuartional_forces/Version_003/__doku_CF3D/cf3d/src/ConF3D/CT_specimen_linear_elastic/CT_specimen_CPE4.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_01 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xy2 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_02 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c2 = session.Curve(xyData=xy2)
xy3 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_03 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c3 = session.Curve(xyData=xy3)
xy4 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_04 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c4 = session.Curve(xyData=xy4)
xy5 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_05 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c5 = session.Curve(xyData=xy5)
xy6 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_06 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c6 = session.Curve(xyData=xy6)
xy7 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_07 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c7 = session.Curve(xyData=xy7)
xy8 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_08 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c8 = session.Curve(xyData=xy8)
xy9 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_09 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c9 = session.Curve(xyData=xy9)
xy10 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_10 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c10 = session.Curve(xyData=xy10)
xy11 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_11 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c11 = session.Curve(xyData=xy11)
xy12 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_12 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c12 = session.Curve(xyData=xy12)
xy13 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_13 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c13 = session.Curve(xyData=xy13)
xy14 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_14 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c14 = session.Curve(xyData=xy14)
xy15 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_15 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c15 = session.Curve(xyData=xy15)
xy16 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_16 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c16 = session.Curve(xyData=xy16)
xy17 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_17 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c17 = session.Curve(xyData=xy17)
xy18 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_18 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c18 = session.Curve(xyData=xy18)
xy19 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_19 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c19 = session.Curve(xyData=xy19)
xy20 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_20 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c20 = session.Curve(xyData=xy20)
xy21 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_21 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c21 = session.Curve(xyData=xy21)
xy22 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_22 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c22 = session.Curve(xyData=xy22)
xy23 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_23 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c23 = session.Curve(xyData=xy23)
xy24 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_24 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c24 = session.Curve(xyData=xy24)
xy25 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_25 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c25 = session.Curve(xyData=xy25)
xy26 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_26 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c26 = session.Curve(xyData=xy26)
xy27 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_27 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c27 = session.Curve(xyData=xy27)
xy28 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_28 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c28 = session.Curve(xyData=xy28)
xy29 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_29 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c29 = session.Curve(xyData=xy29)
xy30 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='J-integral: J at J_ABQ_J_ABQ__PICKEDSET28_Contour_30 in ELSET  ALL ELEMENTS', 
    steps=('Loading', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c30 = session.Curve(xyData=xy30)
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, 
    c12, c13, c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, 
    c27, c28, c29, c30, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
