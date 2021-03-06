# Do not edit this file or it may not load correctly
# if you try to open it with the RSG Dialog Builder.

# Note: thisDir is defined by the Activator class when
#       this file gets exec'd

from rsg.rsgGui import *
from abaqusConstants import INTEGER, FLOAT
dialogBox = RsgDialog(title='Fracture Tool', kernelModule='EDC2', kernelFunction='fEDC', includeApplyBtn=False, includeSeparator=True, okBtnText='OK', applyBtnText='Apply', execDir=thisDir)
RsgTabBook(name='TabBook_1', p='DialogBox', layout='0')
RsgTabItem(name='TabItem_5', p='TabBook_1', text='Instructions')
RsgLabel(p='TabItem_5', text='Before running the program make sure you have:', useBoldFont=False)
RsgLabel(p='TabItem_5', text='- A part model and its mesh (Hex mesh only)', useBoldFont=False)
RsgLabel(p='TabItem_5', text='- Composite definition and its datum axis', useBoldFont=False)
RsgLabel(p='TabItem_5', text='- In the composite definition make sure the relative element thicknesses for the plies add to one', useBoldFont=False)
RsgLabel(p='TabItem_5', text='- In the composite definition make sure only one region is defined for all plies', useBoldFont=False)
RsgLabel(p='TabItem_5', text='-Tool assumes direction one of the composite layup definition is the firbe direction', useBoldFont=False)
RsgLabel(p='TabItem_5', text='- In the stack direction definition, element direction 1 = x-axis, element direction 2 = y-axis, and element direction 3 = z-direction', useBoldFont=False)
RsgLabel(p='TabItem_5', text='Step by step instructions:', useBoldFont=False)
RsgLabel(p='TabItem_5', text='For step by step instructions please refer to the user guide provided with the plugin', useBoldFont=False)
RsgLabel(p='TabItem_5', text='or use the link below for a video demonstration of the tool:', useBoldFont=False)
RsgLabel(p='TabItem_5', text='https://youtu.be/zVOHwa3-SCU', useBoldFont=False)
RsgTabItem(name='TabItem_3', p='TabBook_1', text='Part Selection')
RsgComboBox(name='ComboBox_2', p='TabItem_3', text='Part:', keyword='part', default='', comboType='MDB', repository='parts', rootText='Model:', rootKeyword='model', layout='HORIZONTAL')
RsgTabItem(name='TabItem_4', p='TabBook_1', text='Material Properties')
RsgLabel(p='TabItem_4', text='Define Material', useBoldFont=True)
RsgRadioButton(p='TabItem_4', text='Material from database', keyword='cusMat', default=True)
RsgComboBox(name='ComboBox_3', p='TabItem_4', text='Material Database:', keyword='mat', default='T800s/M21', comboType='STANDARD', repository='', rootText='', rootKeyword=None, layout='')
RsgListItem(p='ComboBox_3', text='T800s/M21')
RsgListItem(p='ComboBox_3', text='T300/920')
RsgListItem(p='ComboBox_3', text='T300/913')
RsgComboBox(name='ComboBox_4', p='TabItem_4', text='Length Unit:', keyword='unit', default='meter', comboType='STANDARD', repository='', rootText='', rootKeyword=None, layout='')
RsgListItem(p='ComboBox_4', text='meter')
RsgListItem(p='ComboBox_4', text='millimeter')
RsgGroupBox(name='GroupBox_1', p='TabItem_4', text='T800s/M21', layout='0')
RsgLabel(p='GroupBox_1', text='Translaminar Fracture Toughness: 209 kJ/m^2', useBoldFont=False)
RsgLabel(p='GroupBox_1', text='Intralaminar Fracture Toughness: 0.255 kJ/m^2', useBoldFont=False)
RsgLabel(p='GroupBox_1', text='Longitudinal Tensile Strength: 3066.96MPa', useBoldFont=False)
RsgGroupBox(name='GroupBox_4', p='TabItem_4', text='T300/920', layout='0')
RsgLabel(p='GroupBox_4', text='Translaminar Fracture Toughness: 133 kJ/m^2', useBoldFont=False)
RsgLabel(p='GroupBox_4', text='Intralaminar Fracture Toughness: 0.456 kJ/m^2', useBoldFont=False)
RsgLabel(p='GroupBox_4', text='Longitudinal Tensile Strength: 3530MPa', useBoldFont=False)
RsgGroupBox(name='GroupBox_5', p='TabItem_4', text='T300/913', layout='0')
RsgLabel(p='GroupBox_5', text='Translaminar Fracture Toughness: 133 kJ/m^2', useBoldFont=False)
RsgLabel(p='GroupBox_5', text='Intralaminar Fracture Toughness: 0.211 kJ/m^2', useBoldFont=False)
RsgLabel(p='GroupBox_5', text='Longitudinal Tensile Strength: 3530MPa', useBoldFont=False)
RsgHorizontalFrame(name='HFrame_3', p='TabItem_4', layout='0', pl=0, pr=0, pt=15, pb=0)
RsgRadioButton(p='TabItem_4', text='Custom Material', keyword='cusMat', default=False)
RsgTextField(p='TabItem_4', fieldType='Float', ncols=12, labelText='Translaminar fracture toughness:', keyword='g0', default='0')
RsgTextField(p='TabItem_4', fieldType='Float', ncols=12, labelText='Intralaminar Fracture Toughness', keyword='g90', default='0')
RsgTextField(p='TabItem_4', fieldType='Float', ncols=12, labelText='Longitudinal Tensile Strength:', keyword='cusSigma', default='0')
RsgTextField(p='TabItem_4', fieldType='Float', ncols=12, labelText='Transverse Tensile Strength:', keyword='cusSigma2', default='0')
RsgTextField(p='TabItem_4', fieldType='Float', ncols=12, labelText='In Plane Shear Strength:', keyword='cusTau', default='0')
RsgLabel(p='TabItem_4', text='Note: Please make sure the defined material properties has units defined consistent witht the FE model', useBoldFont=False)
RsgTabItem(name='TabItem_1', p='TabBook_1', text='Selected Potential Fracture Plane')
RsgCheckButton(p='TabItem_1', text='Estimate energy and strength for desired possible fracture plane', keyword='singlePlane', default=False)
RsgLabel(p='TabItem_1', text='Define desired fracture plane by defining 3 points', useBoldFont=True)
RsgRadioButton(p='TabItem_1', text='Select points in viewport:', keyword='selectPsingle', default=True)
RsgPickButton(p='TabItem_1', text='Point 1', keyword='coord1', prompt='Pick a node', entitiesToPick='MDB_MESH|NODES', numberToPick='ONE')
RsgPickButton(p='TabItem_1', text='Point 2', keyword='coord2', prompt='Pick a node', entitiesToPick='MDB_MESH|NODES', numberToPick='ONE')
RsgPickButton(p='TabItem_1', text='Point 3', keyword='coord3', prompt='Pick a node', entitiesToPick='MDB_MESH|NODES', numberToPick='ONE')
RsgHorizontalFrame(name='HFrame_4', p='TabItem_1', layout='0', pl=0, pr=0, pt=15, pb=0)
RsgRadioButton(p='TabItem_1', text='Define Points Manually:', keyword='selectPsingle', default=False)
RsgVerticalFrame(name='VFrame_2', p='TabItem_1', layout='0', pl=0, pr=0, pt=0, pb=0)
RsgTable(p='VFrame_2', numRows=1, columnData=[('X', 'Float', 100), ('Y', 'Float', 100), ('Z', 'Float', 100)], showRowNumbers=False, showGrids=True, keyword='enterCoord1', popupFlags='')
RsgTable(p='VFrame_2', numRows=1, columnData=[('X', 'Float', 100), ('Y', 'Float', 100), ('Z', 'Float', 100)], showRowNumbers=False, showGrids=True, keyword='enterCoord2', popupFlags='')
RsgTable(p='VFrame_2', numRows=1, columnData=[('X', 'Float', 100), ('Y', 'Float', 100), ('Z', 'Float', 100)], showRowNumbers=False, showGrids=True, keyword='enterCoord3', popupFlags='')
RsgIcon(p='TabItem_1', fileName=r'single plane diagram v2.PNG')
RsgTabItem(name='TabItem_2', p='TabBook_1', text='Sweep Potential Fracture Planes')
RsgCheckButton(p='TabItem_2', text='Sweep Through Potential Fracture Planes', keyword='multiPlane', default=False)
RsgLabel(p='TabItem_2', text='Define sweep axis', useBoldFont=True)
RsgHorizontalFrame(name='HFrame_1', p='TabItem_2', layout='0', pl=0, pr=0, pt=0, pb=0)
RsgRadioButton(p='HFrame_1', text='Use global axis', keyword='selectAxis', default=True)
RsgGroupBox(name='GroupBox_6', p='TabItem_2', text='Global Axis', layout='0')
RsgRadioButton(p='GroupBox_6', text='X', keyword='axisGlobal', default=True)
RsgRadioButton(p='GroupBox_6', text='Y', keyword='axisGlobal', default=False)
RsgRadioButton(p='GroupBox_6', text='Z', keyword='axisGlobal', default=False)
RsgHorizontalFrame(name='HFrame_5', p='TabItem_2', layout='0', pl=0, pr=0, pt=10, pb=0)
RsgRadioButton(p='TabItem_2', text='Use two points', keyword='selectAxis', default=False)
RsgHorizontalFrame(name='HFrame_2', p='TabItem_2', layout='0', pl=0, pr=0, pt=0, pb=0)
RsgGroupBox(name='GroupBox_7', p='TabItem_2', text='Axis Points Definition', layout='0')
RsgRadioButton(p='GroupBox_7', text='Select points in viewport:', keyword='axisPoints', default=True)
RsgPickButton(p='GroupBox_7', text='Start Point', keyword='axisPointStart', prompt='Pick an entity', entitiesToPick='MDB_MESH|NODES', numberToPick='ONE')
RsgPickButton(p='GroupBox_7', text='End Point', keyword='axisPointEnd', prompt='Pick an entity', entitiesToPick='MDB_MESH|NODES', numberToPick='ONE')
RsgRadioButton(p='GroupBox_7', text='Define Points Manually', keyword='axisPoints', default=False)
RsgLabel(p='GroupBox_7', text='Start Coordinates', useBoldFont=False)
RsgTable(p='GroupBox_7', numRows=1, columnData=[('X', 'Float', 100), ('Y', 'Float', 100), ('Z', 'Float', 100)], showRowNumbers=False, showGrids=True, keyword='coordAxisStart', popupFlags='')
RsgLabel(p='GroupBox_7', text='End Coordinates', useBoldFont=False)
RsgTable(p='GroupBox_7', numRows=1, columnData=[('X', 'Float', 100), ('Y', 'Float', 100), ('Z', 'Float', 100)], showRowNumbers=False, showGrids=True, keyword='coordAxisEnd', popupFlags='')
RsgLabel(p='TabItem_2', text='Resolution', useBoldFont=True)
RsgTextField(p='TabItem_2', fieldType='Integer', ncols=12, labelText='Number of potential fracture Planes (N):', keyword='N', default='10')
RsgLabel(p='TabItem_2', text='Note: Planes Will Be Created In The Direction Of Axis unnless specified as in the axis points definition', useBoldFont=False)
RsgIcon(p='TabItem_2', fileName=r'multiplane diagramv2.PNG')
dialogBox.show()