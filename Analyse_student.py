# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
#
# Ceci est un code d'Analyse des resultats obtenues avec le code 'ExerciceY_XXXX'.
# Le code d'analyse est base sur le standard python3 (on recommende python3.6 ou superieur).
# Pour installer les librairies on recommende d'utiliser pip3. Si pip3 n'est pas installe,
# il est possible de suivre la procedure d'un des links suivants:
#   https://linuxconfig.org/how-to-install-pip-on-ubuntu-18-04-bionic-beaver
#   https://linuxize.com/post/how-to-install-pip-on-ubuntu-18.04/
# 
# Ensuite, il faut installer les librairies: 
#   numpy 
#   matplotlib
#   scipy 
#   os
# methode d'installation conseille: utiliser la ligne de commande: 
#   pip3 install --user *nome-librairie*
# dans un terminal linux
#
# Pour utiliser le code d'analyse, il faut que son source (ce ficher) soit 
# dans le repertoire contenant le binaire 'ExerciceY_XXXX' et les 
# fichers d'output '.out'. Pour l'executer, il faut utiliser la ligne de commande
# suivantes dans le terminal linux:
#   python3 Analyse.py
#
# ----------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------
# Modules --------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
import numpy as np
import PlotResults # importer le module pour plotter les resultats
import SchemeAnalysisSuite # importer le module pour analyser les schemas
import ResultAnalysisSuit  # importer le module pour analyser les resultats
import AnalyticalSolutions as ansol # importer le module des solutions analytiques

# ----------------------------------------------------------------------------------------------
# Parametres -----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

# Parametres pour lire et plotter les resultats d'une simulation -------------------------------

# Nome du fichier d'output a analyser
doFigures = False # si vrai les figures pour un fichier fileName sont produites 
delimiter = ' '
fileNames =['']

# dictionaire contenant les parametres pour faire une figure
figure1 = {'NRows':1,'NCols':2} # plot la dynamique
figure1['NPlots'] = [1,1]
figure1['indexes'] = [[[,]],[[,]]]
figure1['titles'] = ['','']
figure1['xLables'] = ['','']
figure1['yLables'] = ['','']
figure1['normalisedPlot'] = [False,False]
figure1['axis'] = ['auto','auto']
figure1['legend'] = [[],[]]
figure1['legendLocation'] = [1,1]
figure1['AnalyticalSolutions'] = []
#figure1['AnalyticalSolutions'] = [ansol.positionPendule,ansol.vitessePendule]

# dictionaire contenant les parametres pour faire une figure
figure2 = {'NRows':1,'NCols':1} # plot la dynamique
figure2['NPlots'] = [1]
figure2['indexes'] = [[[,]]]
figure2['titles'] = ['']
figure2['xLables'] = ['']
figure2['yLables'] = ['']
figure2['normalisedPlot'] = [False]
figure2['axis'] = ['auto']
figure2['legend'] = [[]]
figure2['legendLocation'] = [1]
figure2['AnalyticalSolutions'] = []

# dictionaire contenant les parametres pour faire une figure
figure3 = {'NRows':1,'NCols':1} # plot la dynamique
figure3['NPlots'] = [1]
figure3['indexes'] = [[[,]]]
figure3['titles'] = ['']
figure3['xLables'] = ['']
figure3['yLables'] = ['']
figure3['normalisedPlot'] = [False]
figure3['axis'] = ['auto']
figure3['legend'] = [[]]
figure3['legendLocation'] = [1]
figure3['AnalyticalSolutions'] = []

# dictionaire contenant les parametres pour faire tous figures
figureParameters = {'plotData':[figure1,figure2,figure3]}
figureParameters['fontSize'] = 32
figureParameters['lineStyle'] = ['-','-','-','-','-','-','-','-','-','-','-','-']
figureParameters['lineWidth'] = 3
figureParameters['marker'] = ['.','.','.','.','.','.','.','.','.','.','.','.']
figureParameters['markerSize'] = 12
figureParameters['plotType'] =  '' #'loglog'
#figureParameters['plotType'] =  'loglog'

# Parametres pour lire et plotter en 2D les resultats d'une simulation -------------------------
# Nome du fichier d'output a analyser
do2DSimpleFigures = False # si vrai les figures pour un fichier fileName sont produites 
delimiter = ' '
fileNames2DSimple =['output_psi2.out']

# dictionaire contenant les parametres pour faire une figure
figure1 = {}
figure1['titles'] = 'Fonction onde: oscillateur harmonique'
figure1['xLables'] = 'x'
figure1['yLables'] = 't'

# dictionaire contenant les parametres pour faire tous figures
figure2DSimpleParameters = {'plotData':[figure1]}
figure2DSimpleParameters['fontSize'] = 32
#figure2DSimpleParameters['contour_levels'] = 900 # number of levels contours
figure2DSimpleParameters['contour_levels'] = 300 # number of levels contours
#figure2DSimpleParameters['contour_levels'] = 100 # number of levels contours
#figure2DSimpleParameters['plotType'] = 'contour' # 2D plots contour contourf
figure2DSimpleParameters['plotType'] = 'contourf' # 2D plots contour contourf

# Parametres pour lire et plotter en 2D les resultats d'une simulation -------------------------
# Nome du fichier d'output a analyser
do2DFigures = False # si vrai les figures pour un fichier fileName sont produites 
delimiter = ' '
fileNames2D =['solutionOutput.out']

# dictionaire contenant les parametres pour faire une figure
figure1 = {'NRows':1,'NCols':2} # plot la dynamique
figure1['extractValue'] = [0,591] # values to be extracted
figure1['indexes'] = [0,1,2,3]
figure1['tableSize'] = (81,81)
figure1['titles'] = 'onde'
figure1['xLables'] = 'x'
figure1['yLables'] = 't'

# dictionaire contenant les parametres pour faire tous figures
figure2DParameters = {'plotData':[figure1]}
figure2DParameters['fontSize'] = 32
#figure2DParameters['contour_levels'] = 900 # number of levels contours
figure2DParameters['contour_levels'] = 100 # number of levels contours
#figure2DParameters['plotType'] = 'contour' # 2D plots contour contourf
figure2DParameters['plotType'] = 'contourf' # 2D plots contour contourf

# Parametres pour lire et plotter en 2D les resultats d'une simulation -------------------------
# Nome du fichier d'output a analyser
doVectorFigures = False # si vrai les figures pour un fichier fileName sont produites 
delimiter = ' '
fileNamesVector2D =['fluxOutput.out']

# dictionaire contenant les parametres pour faire une figure
figure1 = {'NRows':1,'NCols':2} # plot la dynamique
figure1['extractValue'] = [0,591] # values to be extracted
# 0:index temps, 1:x, 2:y, 3:flux_x, 4:flux_y, 5:intensite du flux
figure1['indexes'] = [0,1,2,3,4,5,6]
figure1['tableSize'] = (80,80)
figure1['titles'] = 'Density'
figure1['xLables'] = 'x'
figure1['yLables'] = 'y'

# dictionaire contenant les parametres pour faire tous figures
figureVectorParameters = {'plotData':[figure1]}
figureVectorParameters['fontSize'] = 32
figureVectorParameters['plotType'] = ' '

# Parametres pour executer le test de convergence ----------------------------------------------

# Variables definies par l'utilisateur
doConvergenceTest = False # si vrai le test de convergence est execute
# Parametres du fichiers d'input
inputs = {'tfin':0.,'mass':0.,'q':0.,\
'L':0.,'g':0.,'k':0.,'E0':0.,'omega':0.,\
'phi':0.,'theta0':0.,'thetaDot0':0.,'dt':0.,\
'schema':'SV','output':'output.out','sampling':0}
# Parametre pour le test de convergence
convergenceParametres = dict()
convergenceParametres['fileDelimiter'] = ' '
convergenceParametres['typeFigureOfMerit'] = ''
convergenceParametres['arrayType'] = ''
convergenceParametres['selectColumn'] = False
convergenceParametres['doSimulations'] = False
convergenceParametres['binaryFileName'] = ''
convergenceParametres['inputFileName'] = 'configuration.in'
convergenceParametres['doLinearRegression'] = True
convergenceParametres['numberOfPointForRegression'] = 3
convergenceParametres['meshGenerator'] = 'power10'
convergenceParametres['meshValues'] = [0,0,0]
convergenceParametres['meshType'] = np.int32
convergenceParametres['keyForSimulationSeries'] = ''
convergenceParametres['keyForOutputFileName'] = 'output'
convergenceParametres['interpolationValues'] = []
convergenceParametres['interpolationIndexes'] = []
convergenceParametres['indexList'] = []
convergenceParametres['plot1overMesh'] = False
convergenceParametres['analyticalSolution']=[]
convergenceParametres['inputParameters'] = inputs

# parametres pour generer la figure avec l'etude de convergence
figure1 = {'NRows':1,'NCols':2} # plot la dynamique
figure1['NPlots'] = [1,1]
figure1['indexes'] = [[[,]],[[,]]]
figure1['titles'] = ['','']
figure1['xLables'] = ['','']
figure1['yLables'] = ['','']
figure1['normalisedPlot'] = [False,False]
figure1['axis'] = ['auto','auto']
figure1['legend'] = [['slope='],['slope=']]
figure1['legendLocation'] = [1,1]
figure1['AnalyticalSolutions'] = []

#figure2 = {'NRows':1,'NCols':1} # plot la dynamique
#figure2['NPlots'] = [1]
#figure2['indexes'] = [[[,]]]
#figure2['titles'] = ['']
#figure2['xLables'] = ['']
#figure2['yLables'] = ['']
#figure2['normalisedPlot'] = [False]
#figure2['axis'] = ['']
#figure2['legend'] = [['slope=']]
#figure2['legendLocation'] = [1]
#figure2['AnalyticalSolutions'] = []

# dictionaire contenant les parametres pour faire tous figures
convergenceFigureParameters = {'plotData':[figure1]}
convergenceFigureParameters['fontSize'] = 32
convergenceFigureParameters['lineStyle'] = ['-','-','-','-']
convergenceFigureParameters['lineWidth'] = 3
convergenceFigureParameters['marker'] = ['x','x','x','x']
convergenceFigureParameters['markerSize'] = 18
convergenceFigureParameters['plotType'] = 'loglog'

# Parametres pour executer etudes parametriques ----------------------------------------------
# Variables definies par l'utilisateur
doParametricStudy = False # si vrai le test de convergence est execute
# Parametres du fichiers d'input
inputs = {'tfin':0.,'mass':0.,'q':0.,\
'L':0.,'g':0.,'k':0.,'E0':0.,'omega':0.,\
'phi':0.,'theta0':0.,'thetaDot0':0.,'dt':0.,\
'schema':'SV','output':'output.out','sampling':0}
# Parametre pour l'etude parametrique
studyParametres = dict()
studyParametres['fileDelimiter'] = ' '
studyParametres['doSimulations'] = True
studyParametres['binaryFileName'] = ''
studyParametres['inputFileName'] = 'configuration.in'
studyParametres['meshGenerator'] = ''
studyParametres['meshValues'] = [0.,0.,0.] 
studyParametres['meshType'] = np.int32
studyParametres['keyForSimulationSeries'] = ''
studyParametres['keyForOutputFileName'] = 'output'
studyParametres['action'] = ['']
studyParametres['indexList'] = [3]
studyParametres['inputParameters'] = inputs

# parametres pour generer la figure avec l'etude de convergence
figure1 = {'NRows':1,'NCols':1} # plot la dynamique
figure1['NPlots'] = [1]
figure1['indexes'] = [[[,]]]
figure1['titles'] = ['']
figure1['xLables'] = ['']
figure1['yLables'] = ['']
figure1['normalisedPlot'] = [False]
figure1['axis'] = ['auto']
figure1['legend'] = [[]]
figure1['legendLocation'] = [1]
figure1['AnalyticalSolutions'] = []

# dictionaire contenant les parametres pour faire tous figures
studyFigureParameters = {'plotData':[figure1]}
studyFigureParameters['fontSize'] = 32
studyFigureParameters['lineStyle'] = ['-']
studyFigureParameters['lineWidth'] = 3
studyFigureParameters['marker'] = ['x']
studyFigureParameters['markerSize'] = 12
studyFigureParameters['plotType'] = ' '

# Parametres pour faire des analyse des donnees ------------------------------------------------
doResultAnalysis = False # executer les analyses si True
analysisFileNames = ['']
analysisDelimiter = ' '

# parametres pour analyser les distances
distanceParameters = {} # initialiser le dictionaire
# sauvegarder l'index de la mesh
distanceParameters['meshIndex'] = 0
# sauvegarder la list des indexes pour chaque distance
distanceParameters['indexes'] = [[,]]
# sauvegarder la list des normalisations pour chaque distance
distanceParameters['normalisation'] = [[]]
# dictionaire avec les analyse a faire
analysisParameters = {'distance':distanceParameters}

# parametres pour generer la figure avec l'etude de convergence
figure1 = {'NRows':1,'NCols':1} # plot la dynamique
figure1['NPlots'] = []
figure1['indexes'] = [[[]]]
figure1['titles'] = ['']
figure1['xLables'] = ['']
figure1['yLables'] = ['']
figure1['normalisedPlot'] = [False]
figure1['axis'] = ['auto']
figure1['legend'] = [[]]
figure1['legendLocation'] = [1]
figure1['AnalyticalSolutions'] = []

# comparaison puissance forces NC et mecanique
#puissanceParameters = {}
#puissanceParameters['derivativeType'] = 'forwardDerivative' # type de derivee
#puissanceParameters['meshIndex'] = 0 # index du maillage
#puissanceParameters['index'] = 3 # index de l'energie mecanique
#puissanceParameters['comparisonIndex'] = 4 # index de la puissance NC
#analysisParameters = {'computeCompareDerivative1D':puissanceParameters}

#figure1 = {'NRows':1,'NCols':1} # plot la dynamique
#figure1['NPlots'] = []
#figure1['indexes'] = [[[]]]
#figure1['titles'] = ['']
#figure1['xLables'] = ['']
#figure1['yLables'] = ['']
#figure1['normalisedPlot'] = [False]
#figure1['axis'] = ['auto']
#figure1['legend'] = [['','']]
#figure1['legendLocation'] = [1]
#figure1['AnalyticalSolutions'] = []


# dictionaire contenant les parametres pour faire tous figures
analysisFigureParameters = {'plotData':[figure1]}
analysisFigureParameters['fontSize'] = 32
analysisFigureParameters['lineStyle'] = ['-','-','-','-']
analysisFigureParameters['lineWidth'] = 3
analysisFigureParameters['marker'] = ['x','x','x','x']
analysisFigureParameters['markerSize'] = 0
analysisFigureParameters['plotType'] = 'semilogy'

# ----------------------------------------------------------------------------------------------
# Lire et plotter les resultats de la simulation -----------------------------------------------
# ----------------------------------------------------------------------------------------------

# verifier s'il faut faire les figures
if(doFigures):
  # Chargement des donnees dans un numpy array
  # initialiser la class pour faire des plots
  plotResult = PlotResults.PlotResults(fileNameList=fileNames,figureParameters=figureParameters,\
  valuesDelimiter=delimiter)
  # faire des plots simples
  plotResult.SimplePlotFigures()

# verifier s'il faut faire les figures 2D en utilisant des tableu
if(do2DSimpleFigures):
  # Chargement des donnees dans un numpy array
  # initialiser la class pour faire des plots
  plotResult = PlotResults.PlotResults(fileNameList=fileNames2DSimple,\
  figureParameters=figure2DSimpleParameters,\
  valuesDelimiter=delimiter)
  # faire des plots simples
  plotResult.plotSimple2D()

# verifier s'il faut faire les figures 2D
if(do2DFigures):
  # Chargement des donnees dans un numpy array
  # initialiser la class pour faire des plots
  plotResult = PlotResults.PlotResults(fileNameList=fileNames2D,\
  figureParameters=figure2DParameters,\
  valuesDelimiter=delimiter)
  # faire des plots simples
  plotResult.plotSimple2DTimeFlattenedFigure()

# verifier s'il faut faire les figures avec vecteurs
if(doVectorFigures):
  # Chargement des donnees dans un numpy array
  # initialiser la class pour faire des plots
  plotResult = PlotResults.PlotResults(fileNameList=fileNamesVector2D,\
  figureParameters=figureVectorParameters,\
  valuesDelimiter=delimiter)
  # faire des plots simples avec vecteurs
  plotResult.plotSimple2DTimeFlattenedVectorFigure()


# ----------------------------------------------------------------------------------------------
# Executer le test de convergence --------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

# verifier qu'il faut faire le test de convergence
if(doConvergenceTest):

  # initialise un objet SchemeAnalysisSuite
  schemeAnalysis = SchemeAnalysisSuite.SchemeAnalysisSuite(\
  convergenceParameters=convergenceParametres,\
  figureParameters=convergenceFigureParameters)

  # executer le test de convergence
  schemeAnalysis.doConvergenceTest()

# ----------------------------------------------------------------------------------------------
# Executer les etudes parametriques ------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

# verifier qu'il faut faire les etudes parametriques
if(doParametricStudy):

  # initialise un objet SchemeAnalysisSuite
  schemeAnalysis = SchemeAnalysisSuite.SchemeAnalysisSuite(\
  convergenceParameters=studyParametres,\
  figureParameters=studyFigureParameters)

  # executer les etudes parametriques
  schemeAnalysis.doParameterStudy()

# ----------------------------------------------------------------------------------------------
# Executer les analyses des resultats ----------------------------------------------------------
# ----------------------------------------------------------------------------------------------

# verifier qu'il faut faire des analyses des resultats
if(doResultAnalysis):

  # initialise un objet ResultsAnalysisSuite
  resultAnalysis = ResultAnalysisSuit.ResultAnalysisSuit(fileNameList=analysisFileNames,\
                   parameters=analysisParameters,\
                   figureParameters=analysisFigureParameters,\
                   valuesDelimiter=analysisDelimiter)

  # executer les analyses
  resultAnalysis.doResultAnalysis()

# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

