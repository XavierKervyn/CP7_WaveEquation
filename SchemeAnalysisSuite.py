# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
#
# Ceci est un module du code d'Analyse des schemas numeriques qui a comme objectifs
# d'etudier leur convergence et proprietes de conservation.
# Le code d'analyse est base sur le standard python3 (on recommende python3.6 ou superieur).
# Pour installer les librairies on recommende d'utiliser pip3. Si pip3 n'est pas installe,
# il est possible de suivre la procedure d'un des links suivants:
#   https://linuxconfig.org/how-to-install-pip-on-ubuntu-18-04-bionic-beaver
#   https://linuxize.com/post/how-to-install-pip-on-ubuntu-18.04/
# 
# Ensuite, il faut installer les librairies: 
#   numpy 
#   scipy
#   os
# methode d'installation conseille: utiliser la ligne de commande: 
#   pip3 install --user *nome-librairie*
# dans un terminal linux
#
# Pour utiliser le code d'analyse, il faut que son source (ce ficher) soit 
# dans le repertoire contenant le main Analyse.py. 
#
# ----------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------
# Librairies -----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
import numpy as np # importer numpy pour tout le module

# creer la class SchemeAnalysisSuite
class SchemeAnalysisSuite():

  # constructeur de la class
  def __init__(self,convergenceParameters={},figureParameters={}):

    # sauvegarder les parametres de convergence
    # Ce dictionaire contient tout les parametres pour
    # effectuer un etude de convergence ou conservation
    # parametres:
    #   typeFigureOfMerit:	    (str) for array2D select the figure of merit to
    #					  compute. Available:
    #					    array1DAbs: compute the absolute of the error
    #					    for 1D array 
    #					    integral1D: compute a 1D integral from
    #					    of the error 2D arrays via trapezoidal integral
    #   arrayType:                  (str)  type of array:
    #					   array1D: array with only one solution dimension
    #					   array2D: array with two solution dimensions but
    #					   with only one solution
    #					   flattenedArray2D: array with two or more
    #					   solution dimension and with one or more
    #					   solutions flattened in a 2D array. The
    #					   column order must be: dimension-N-1,
    #					   dimension-N-2,...,dimension-1,dimension-0,
    #					   solution-0,solution-1,...,solution-N
    #   selectColumn:	            (bool) if true figure of merits are computed
    #					   column-wise
    #   fileDelimeter:		    (str) character qui delimite les colonnes
    #					  du fichier d'output
    #   binaryFileName:		    (str) nom de l'executable
    #   inputFileName:		    (str) nom du fichier d'input
    #   numberOfPointForRegression: (int) nombre de points pour effecture
    #					  une regressione lineaire
    #   meshGenerator:		    (str) nombre de points pour generer un
    #					  maillage, available:
    #					  0) input
    #					  1) linear
    #					  2) power10
    #   meshValues:		    (list(float64))(NValues) valeurs pour
    #					  generer le maillage: si input
    #					  est utilise meshValues doit 
    #					  contenir le maillage entier
    #					  autrement seulement les valeurs
    #					  minimale et maximale sont demandees
    #					  le nombre de points est facultatif
    #   meshType:		    (dataType) mesh dtype
    #   keyForSimulationSeries:	    (str) clee du parametre de convergence
    #					  (nom du parametre utilise pour
    #					  l'etude de convergence)
    #   keyForOutputFileName:	    (str) clee du parametre contenant le nom
    #					  du fichier d'output
    #   interpolationValues:	    (list(float64)) abscisse(s) d'interpolation
    #					  pour tous les valeurs
    #   interpolationIndexes:	    (list(int)) indexe(s) des abscisse(s)
    #					  d'interpolation pour tous les valeurs
    #   indexList:		    (list(int))(NIndexes) liste des indexes des
    #					  ordonees utilisees a etudier
    #   analyticalSolution:	    (list(float64)) liste des solutions analytique
    #					  avec lesquelles calculer l'erreur 
    #					  de convergence ou des fonctions
    #					  calculant les solution analytique
    #   inputParameters:	    (dict) dictionaire avec les parametres
    #					   a ecrire dans le ficher d'input
    #   action:           	    (dict) dictionaire qui definisse le type
    #					   d action a faire pour un etude (cle)
    #					   et ses parametres (0 si aucun)
    #   doSimulations:		    (bool) if True new simulations are run
    self.convergenceParameters = convergenceParameters

    # sauvegarder les parametres pour faire les figures
    # attention: la methode de doConvergenceTest demande
    # un figure avec NIndexes subplots. Les parametres
    # de la figure de convergence doit etre le premier
    # dictionaire dans figureParameters['plotData']
    # pour plus d'information: regarder le module PlotResults
    self.figureParameters = figureParameters

    # liste des nomes du ficher d'output
    self._outputFileNameList = []
    # array de valeurs utilisees pour etudier les schema
    # comme par, exemple, les valeur de pas de temps pour
    # l'etude de convergence. Taille: (NMesh) 
    self.__mesh = np.array([])
    # array des resultats. Taille: (NIndexes,NMesh)
    self.__results = np.array([])
    # array des erreurs
    self.__errors = np.array([])
    # array des valuers pour les etudes parametrques
    self.__studyValues = np.array([])
    # array des regressions lineaires des erreurs
    # 0) slopes 1) intersection
    self.__errorsLinearRegression = \
    np.array((2,len(self.convergenceParameters['indexList'])),dtype=np.float64)
    # liste contenant les coordonnees des nodes d'interpolation
    # [[coordonnees node 1],[coordonnees node 2],...]
    self.__interpolationNodeCoordinates = []
    # liste contenant les indexes des nodes d'interpolation
    # [[indexes node 1],[indexes node 2],...]
    self.__interpolationNodeIndexes = []


# ----------------------------------------------------------------------------------------------

  # destructeur de la class
  def __del__(self):

      # imprimer le nome de la class
      print(self.__class__.__name__ ,'is deleted')

# ----------------------------------------------------------------------------------------------

  # cette methode retourne la mesh
  # output:
  #   mesh: (array)(NMesh) array contenant le maillage
  def getMesh(self):
  
    # retourner le maillage
    return self.__mesh

# ----------------------------------------------------------------------------------------------

  # cette methode permets de changer le maillage
  # inputs:
  #   mesh: (array)(NMesh) array contenant le maillage
  def setMesh(self,mesh):

    # copier le maillage
    self.__mesh = np.copy(mesh)

# ----------------------------------------------------------------------------------------------

  # cette methode retourne la resultats
  # output:
  #   results: (array)(NMesh,NIndexes) array contenant les resultats
  def getResults(self):
  
    # retourner le maillage
    return self.__results

# ----------------------------------------------------------------------------------------------

  # Cette methode construit une maillage
  def generateMesh(self):

    # verifier la procedure pour genere une mesh
    if(self.convergenceParameters['meshGenerator']=='input'):
      # copier le maillage
      self.__mesh = np.copy(self.convergenceParameters['meshValues'])
    else:
      # verifier si le nombre de points est donne
      if(len(self.convergenceParameters['meshValues'])==2):
        # ajuter le nombre de points
        self.convergenceParameters['meshValues'].append(\
        self.convergenceParameters['meshValues'][1]-\
        self.convergenceParameters['meshValues'][0]+1)
      # calculer le maillage lineaire
      self.__mesh = np.linspace(self.convergenceParameters['meshValues'][0],\
      self.convergenceParameters['meshValues'][1],\
      num=self.convergenceParameters['meshValues'][2],\
      endpoint=True,dtype=np.float64)
      # verifier le type de maillage
      if(self.convergenceParameters['meshGenerator']=='power10'):
        # construire un maillage exposant10
        self.__mesh = np.power(10.0,self.__mesh,dtype=np.float64)
      if(self.convergenceParameters['meshGenerator']=='power2'):
        # construire un maillage exposant10
        self.__mesh = np.power(2.0,self.__mesh,dtype=np.float64)
      if(self.convergenceParameters['meshGenerator']=='square'):
        # construire un maillage exposant10
        self.__mesh = np.power(self.__mesh,2.0,dtype=np.float64)
      if("subtractMeshToValue" in self.convergenceParameters):
        self.__mesh = self.convergenceParameters["subtractMeshToValue"]-self.__mesh
      if("addMeshToValue" in self.convergenceParameters):
        self.__mesh = self.convergenceParameters["addMeshToValue"]+self.__mesh
      # cast mesh datatype
      self.__mesh.astype(self.convergenceParameters['meshType'])

# ----------------------------------------------------------------------------------------------

  # cette methode permet de executer le binaire 
  # si (doSimulations=True) et lire les donnees.
  # inputs:
  #   meshValue: (value) valeur de l'element du maillage
  #			 a ecrire dans le fichier d'input
  def runExecutableAndReadData(self,meshValue):

    # importer la librerie os
    import os
  
    # surecrire la valeur de convergence pour le fichier # d'input
    sign_mesh = self.convergenceParameters['inputParameters'][\
    self.convergenceParameters['keyForSimulationSeries']]/abs(\
    self.convergenceParameters['inputParameters'][\
    self.convergenceParameters['keyForSimulationSeries']]) 
    self.convergenceParameters['inputParameters'][\
    self.convergenceParameters['keyForSimulationSeries']]=sign_mesh*meshValue
    # generer un nome pour le fichier d'output
    self.convergenceParameters['inputParameters'][\
    self.convergenceParameters['keyForOutputFileName']]=\
    ''.join(['output_',\
    self.convergenceParameters['keyForSimulationSeries'],str(meshValue),'.out'])
    # ajouter le nome du fichier d'output a la list
    self._outputFileNameList.append(\
    self.convergenceParameters['inputParameters'][\
    self.convergenceParameters['keyForOutputFileName']])
    # open input file
    with open(self.convergenceParameters['inputFileName'],'w') as inputFile:
      # boucle sur les clees du fichier d'input
      for key,value in self.convergenceParameters['inputParameters'].items():
        # ecrire valeurs dans le fichier d'input
        inputFile.write(''.join([key,' = ',str(value),'\n']))
    # check whether new simulations should be performed
    if(self.convergenceParameters['doSimulations']):
      # executer le binaire
      os.system(''.join(['./',self.convergenceParameters['binaryFileName'],\
      ' ',self.convergenceParameters['inputFileName']]))
    # print
    print('Reading file: ',self.convergenceParameters['inputParameters'][\
    self.convergenceParameters['keyForOutputFileName']])
    # lire les resultats
    self.__results = np.loadtxt(\
    self.convergenceParameters['inputParameters'][\
    self.convergenceParameters['keyForOutputFileName']],\
    dtype=np.float64,delimiter=self.convergenceParameters['fileDelimiter'])  

# ----------------------------------------------------------------------------------------------

  # cette methode selectionne la methode pour chercher les nodes
  # d'interpolation et celle d'interpolation en fonction du nombre 
  # de coordonnees.
  # outputs:
  #   interpolator: (func) methode d'interpolation
  #   nodeFinder:   (func) methode pour trouver les nodes d'interpolation
  def selectInterpolationAndNodeFinderMethods(self):

    interpolator=None # initialiser la methode d'interpolation
    nodeFinder=None   # initialiser la methode pour trouver les nodes
    # verifier la taille de la list des index d'interpolation
    if(len(self.convergenceParameters['interpolationIndexes'])==1):
      # copier la methode d'interpolation
      interpolator = self.interpolator1D
      # copier la methode pour trouver les nodes d'interpolation
      nodeFinder = self.findInterpolationNodes1D
    else:
      # imprimer un message d'erreur
      print('Error: interpolation scheme not implemented')

    # retourner les deux methodes
    return interpolator,nodeFinder

# ----------------------------------------------------------------------------------------------

  # cette methode permet d'interpoler valeurs en 1D
  # inputs: 
  #   valueIndex: (int) index du valeur qu'il faut interpoler
  # outputs:
  #   value: (float64) interpolated value
  def interpolator1D(self,valueIndex):

    # preparer les donnees
    x = [coodinate[0] for coodinate in self.__interpolationNodeCoordinates]
    y = [self.__results[y_value[0],valueIndex] for y_value in self.__interpolationNodeIndexes]
    if(not (self.convergenceParameters['interpolationValues'][0]>=np.amin(x) and\
      self.convergenceParameters['interpolationValues'][0]<=np.amax(x))):
      print("Warning value not in interpolation interval")

    # effectuer l'interpolation
    if(len(self.__interpolationNodeIndexes)>2):
      # utilser les interpolations quadratiques
      from scipy import interpolate as interp
      x = [coodinate[0] for coodinate in self.__interpolationNodeCoordinates]
      y = [self.__results[y_value[0],valueIndex] for y_value in self.__interpolationNodeIndexes]
      f = interp.interp1d(x,y,kind='quadratic',bounds_error=False)
      value = f(self.convergenceParameters['interpolationValues'][0])
    else:
      # utiliser les interpolations lineaires
      value = np.interp(self.convergenceParameters['interpolationValues'][0],x,y)

    return value

# ----------------------------------------------------------------------------------------------

  # cette methode recherche le point du maillage (1D) 
  # dans l'array des solutions et interpole linearment
  # la solution. Si la valeur du maillage n'est pas
  # trouve le dernier resultat est retourne
  def findInterpolationNodes1D(self):

    # chercher l'index le plus proche du valeur donnee
    nearestNodeIndex = (np.abs(\
    self.__results[:,self.convergenceParameters['interpolationIndexes'][0]]-\
    self.convergenceParameters['interpolationValues'][0])).argmin()
    # extraire le valeur du node le plus proche
    nearestInterpolationValue = self.__results[\
    nearestNodeIndex,self.convergenceParameters['interpolationIndexes'][0]]
    # verifier s'il est le node de gauche
    if(nearestInterpolationValue<=self.convergenceParameters['interpolationValues'][0]):
      # copier le nearestNodeIndex comme premier node et nearestNodeIndex+1 comme deuxieme
      if(nearestNodeIndex==(self.__results.shape[0]-1)):
        nearestNodeIndex = nearestNodeIndex-1
      self.__interpolationNodeIndexes=[[nearestNodeIndex],[nearestNodeIndex+1]]
    else:
      # copier le nearestNodeIndex-1 comme premier node et nearestNodeIndex comme deuxieme
      if(nearestNodeIndex==0):
        nearestNodeIndex = 1
      self.__interpolationNodeIndexes=[[nearestNodeIndex-1],[nearestNodeIndex]]
    # ajouter un troisieme point si possible pour interpolation parabolique
    if(self.__interpolationNodeIndexes[1][0]+1 < len(self.__results[:,0])):
       self.__interpolationNodeIndexes.append([self.__interpolationNodeIndexes[1][0]+1])
    elif(self.__interpolationNodeIndexes[0][0]-1 >= 0):
       self.__interpolationNodeIndexes.insert(0,[self.__interpolationNodeIndexes[0][0]-1])
    # boucle sur les valeurs au nodes
    self.__interpolationNodeCoordinates = []
    for nodeCoordinates in self.__interpolationNodeIndexes:
      # sauvegarder les valeurs au nodes
      self.__interpolationNodeCoordinates.append([\
      self.__results[nodeCoordinates[0],self.convergenceParameters['interpolationIndexes'][0]]])

# ----------------------------------------------------------------------------------------------

  # cette methode retourne la fonction de 
  # regressione lineaire desideree en fonction
  # du type du maillage
  # output:
  #   regression: (funct) regression lineaire
  def selectLinearRegression(self):

    # importer la librerie avec la regression lineaire 
    import scipy.stats.mstats as mstats

    # select la regression lineare
    if(self.convergenceParameters['meshGenerator']=='power10'):
      # utiliser la regression du logarithm base 10
      return self.linearRegressionLog10
    else:
      # utiliser la regression lineaire standard
      return mstats.linregress

# ----------------------------------------------------------------------------------------------

  # cette methode overload la regression lineare avec
  # celui du logarithme
  # input:
  #   a: (value) premier valeur
  #   b: (value) deuxieme valeur
  # output:
  #   meme que la refression lineaire
  def linearRegressionLog10(self,a,b):

    # importer la librerie avec la regression lineaire 
    import scipy.stats.mstats as mstats

    # retourner la regression base 10
    return mstats.linregress(np.log10(np.abs(a)),np.log10(np.abs(b)))
  
# ----------------------------------------------------------------------------------------------

  # cette methode execute l'etude de convergence pour une serie de valeurs
  def doConvergenceTest(self):

    # importer le module pour faire des figures
    import PlotResults

    # calculer le maillage pour le test de convergence
    self.generateMesh()
    # check if a flattened 2D array has not to be used
    if(self.convergenceParameters['arrayType']=='array2D'):
      # override number of list
      self.convergenceParameters['indexList'] = [1]
    # initialiser l'array des erreurs
    self.__errors = np.zeros((self.__mesh.shape[0],\
    len(self.convergenceParameters['indexList'])),dtype=np.float64)

    # extract method for computing the figure of merits
    computeFigureOfMerit,interpolator,findNode = self.extractComputeFigureOfMerit()

    # boucle sur les elements du maillage
    for elementId,element in enumerate(self.__mesh):
      # executer une nouvelle simulation
      self.runExecutableAndReadData(element)
      # boucle sur les indexes des valeurs de convergence
      for indexId,index in enumerate(self.convergenceParameters['indexList']):
        # check if a flattened 2D array has not to be used
        if(self.convergenceParameters['arrayType']=='flattenedArray2D'):
          print("method for restoring arrays not implemented yet")
        else:
          # use the result as array2D
          array = self.__results
        # compute figure of merit
        self.__errors[elementId,indexId] = computeFigureOfMerit(indexId,\
        index,array,interpolator,findNode)

    # if linear regression test is required
    if(self.convergenceParameters['doLinearRegression']):
      # select linear regression method
      linearRegression = self.selectLinearRegression()
      # boucle sur les index
      for indexId,index in enumerate(self.convergenceParameters['indexList']):
        # calculer la regression lineaire sur les derniers x-points
        slope,intersection,r_value,p_value,std_err = \
        linearRegression(self.__mesh[self.__mesh.shape[0]-\
        self.convergenceParameters['numberOfPointForRegression']:-1],\
        self.__errors[self.__mesh.shape[0]-\
        self.convergenceParameters['numberOfPointForRegression']:-1,indexId])
        # imprimer la pente
        print('Convergence test index: ',index,' slope: ',slope)
        # verifier la presence de la figure pour la convergence
        if(len(self.figureParameters['plotData'])!=0):
          # sauvegarder la pente comme legende
          self.figureParameters['plotData'][0]['legend'].append(\
          [''.join(['slope: ',str(np.around(slope,decimals=2))])])

    # verifier la presence de la figure pour la convergence
    if(len(self.figureParameters['plotData'])!=0):
      # creer le dictionaire des parametres
      parameters = self.figureParameters.copy()
      # remove all other plot data
      parameters['plotData']=[self.figureParameters['plotData'][0]]
      # initialiser array pour plot
      plotValues = np.zeros((self.__errors.shape[0],\
      self.__errors.shape[1]+1),dtype=np.float64)
      # copier la mesh
      if(not self.convergenceParameters['plot1overMesh']):
        plotValues[:,0] = self.__mesh
      else:
        plotValues[:,0] = 1./self.__mesh 
      # copier les erreurs
      plotValues[:,1:] = self.__errors
      # initialiser un objet pour generer des figures
      plotResult = PlotResults.PlotResults(values=[plotValues],\
      figureParameters=parameters)
      # fare la figure
      plotResult.SimplePlotFigures()

# ----------------------------------------------------------------------------------------------

  # cette methode retourne la methode pour calculer la figures de merite choisit
  def extractComputeFigureOfMerit(self):

    # select the method
    if(self.convergenceParameters['typeFigureOfMerit']=='array1DAbs'):
      # initialiser les methodes d'interpolation
      # et pour pour trouver les nodes d'interpolation
      interpolator,findNode = self.selectInterpolationAndNodeFinderMethods()
      return self.errorAbsArray1D,interpolator,findNode
    elif(self.convergenceParameters['typeFigureOfMerit']=='array1DMaxAbs'):
      return self.errorMaxAbsArray1D,None,None
    elif(self.convergenceParameters['typeFigureOfMerit']=='array2Dintegral1D'):
      return self.intergral1DArray2D,None,None
    else:
      print("Error in SchemeAnalysisSuite: figure of merits not defined")

# ----------------------------------------------------------------------------------------------
  # calculer la norm infinie de l erreur numerique
  # pour des array 1D
  def errorMaxAbsArray1D(self,valueIndexId,valueIndex,array2D,interpolator,findNode):

    # calculer ou extraire la solution analytique
    analyticalSolution = 0.
    if('analyticalSolution' in self.convergenceParameters):
      # verifier que l element soit une methode
      if(callable(self.convergenceParameters['analyticalSolution'][valueIndexId])):
        # calculer la solution analytique
        analyticalSolution = self.convergenceParameters['analyticalSolution'][valueIndexId](\
        array2D[:,self.convergenceParameters['interpolationIndexes'][0]]) 
      else:
        # extraire la valeur de la solution analytique
        analyticalSolution = self.convergenceParameters['analyticalSolution'][valueIndexId]

    # calculer l erreur  
    error = np.amax(np.abs(array2D[:,valueIndex]-analyticalSolution))  
    # returner l erreur
    return error


# ----------------------------------------------------------------------------------------------
  # calculer la valeur absolue de l erreur numerique
  # pour des array 1D
  def errorAbsArray1D(self,valueIndexId,valueIndex,array2D,interpolator,findNode):

    # extraire les nodes d'interpolation
    findNode() 

    # calculer ou extraire la solution analytique
    analyticalSolution = 0.
    if('analyticalSolution' in self.convergenceParameters):
      # verifier que l element soit une methode
      if(callable(self.convergenceParameters['analyticalSolution'][valueIndexId])):
        # calculer la solution analytique
        analyticalSolution = self.convergenceParameters['analyticalSolution'][valueIndexId](\
        np.array(self.convergenceParameters['interpolationValues'][valueIndexId])) 
      else:
        # extraire la valeur de la solution analytique
        analyticalSolution = self.convergenceParameters['analyticalSolution'][valueIndexId]
  
    # calculer l erreur  
    error = np.abs(interpolator(valueIndex)-analyticalSolution)  
    # returner l erreur
    return error

# ----------------------------------------------------------------------------------------------

  # calculer l integrale 1D en utilisant la methode des trapezes
  # par default les tableau sont considere compose comme:
  # valeur [0,0]: valeur de garde
  # premier ligne: valeurs de l abscisse
  # premier colonne: valeurs de l ordonne
  def intergral1DArray2D(self,valueIndexId,valueIndex,array2D,interpolator,findNode):

    # verifier si on doit integrer sur les colonnes ou no
    if(self.convergenceParameters['selectColumn']):
      # en integre sur les colonnes et recherche sur les lignes
      integralAxis = array2D[1:,0]
      researchAxis = array2D[0,1:]
    else:
      # en integre sur les lignes et recherche sur les colonnes
      integralAxis = array2D[0,1:]
      researchAxis = array2D[1:,0]

    # rechercher la valeur plus proches a celle choisit
    integrationValueId = np.argmin(np.abs(researchAxis-
    self.convergenceParameters['interpolationValues'][valueIndexId]))

    # extract integration arrays
    if(self.convergenceParameters['selectColumn']):
      integrationSet = array2D[1:,valueIndexId+1]
    else:
      integrationSet = array2D[valueIndexId+1,1:]

    # verifier la presence d une solution analytique
    if('analyticalSolution' in self.convergenceParameters):
      # check if the analytical solution is a funcion
      if(callable(self.convergenceParameters['analyticalSolution'][valueIndexId])):
        # subtract the analytical solution
        integrationSet = np.abs(integrationSet - \
        self.convergenceParameters['analyticalSolution'][valueIndexId](\
        np.array(integralAxis),\
        np.array(self.convergenceParameters['interpolationValues'][valueIndexId])))
      else:
        # the analytical solution is a list
        analyticalSet = np.array(self.convergenceParameters['analyticalSolution'][valueIndexId])
        # check consistency
        if(np.all(integrationSet.shape==analyticalSet.shape)):
          # compute the error
          integrationSet = np.abs(integrationSet - analyticalSet)
        else:
          print('Error in integral1DArray2D: analytical and numerical solution have different shapes!')
          return 

    # intialise the integration procedure
    integral = 0.e0
    for elementId,element in enumerate(integrationSet[1:]):
      integral = integral + 0.5*(element+integrationSet[elementId])*\
      (integralAxis[elementId+1]-integralAxis[elementId])
    
    # return integral
    return integral

# ----------------------------------------------------------------------------------------------

  # cette methode execute un etude parametrique pour une serie de valeurs
  def doParameterStudy(self):

    # importer le module pour faire des figures
    import PlotResults

    # calculer le maillage pour l etude parametrique
    self.generateMesh()
    # initialiser les actions a faire
    actionsList = self.action()
    # extrare le nombre indexes a faire
    NIndexes = len(self.convergenceParameters['indexList'])
    # initialiser l array contenant les valeurs pour chaque action et index
    self.__studyValues = np.zeros((self.__mesh.shape[0],\
    NIndexes*len(self.convergenceParameters['action'])),dtype=np.float64)

    # creer le dictionaire des parametres
    parameters = self.figureParameters.copy()
    # remove all other plot data
    parameters['plotData']=[self.figureParameters['plotData'][0]]
    # reset the number of lines per subplot
    parameters['plotData'][0]['NPlots'] = []
    # reset the set of indexes
    parameters['plotData'][0]['indexes'] = []

    # boucle sur les elements du maillage
    for elementId,element in enumerate(self.__mesh):
      # executer une nouvelle simulation
      self.runExecutableAndReadData(element)
      # boucle sur les actions
      for actionId,action in enumerate(actionsList):
        # boucle sur les indexes
        for indexId,index in enumerate(self.convergenceParameters['indexList']):
          # executer l'action sur une colonne
          self.__studyValues[elementId,NIndexes*actionId+actionId]=\
          action(self.__results[:,index])
    # boucle sur les actions
    for actionId,action in enumerate(actionsList):
      # add the new number of lines per subplot
      parameters['plotData'][0]['NPlots'].append(NIndexes)
      # initialise a sequence of indexes per subplots
      subplotIndexes = []
      # boucle sur les indexes
      for indexId,index in enumerate(self.convergenceParameters['indexList']):
        # append the new index
        subplotIndexes.append([0,NIndexes*actionId+indexId+1])
        # append the subplot list of indexes
        parameters['plotData'][0]['indexes'].append(subplotIndexes)

    # Plot figures
    # initialiser array pour plot
    plotValues = np.zeros((self.__studyValues.shape[0],\
    self.__studyValues.shape[1]+1),dtype=np.float64)
    # copier la mesh
    plotValues[:,0] = self.__mesh
    # copier les erreurs
    plotValues[:,1:] = self.__studyValues
    # initialiser un objet pour generer des figures
    plotResult = PlotResults.PlotResults(values=[plotValues],\
    figureParameters=parameters)
    # fare la figure
    plotResult.SimplePlotFigures()


# ----------------------------------------------------------------------------------------------

  # cette methode selectionne une action pour effectuer un etude
  def action(self):
    # import modules
    import MathematicalTools as mathTools
    # initialiser la list des actions a faire
    actionsList = [];
    # verifier que la variable action est definite (longueur!=0)
    if(len(self.convergenceParameters['action'])!=0):
      # verifier si le max cest une action a faire
      if('max' in self.convergenceParameters['action']):
        # add numpy np.amax()
        actionsList.append(np.amax)
      # verifier si le max(abs) cest une action a faire
      if('maxabs' in self.convergenceParameters['action']):
        # add numpy np.amax(np.abs())
        actionsList.append(np.amax(np.abs))
      # verifier s il faut chercher le premier pique
      if('first_peak_value' in self.convergenceParameters['action']):
        # utiliser la fonction first_peak
        actionsList.append(mathTools.first_peak_value)
    # return the function list
    return actionsList

# ----------------------------------------------------------------------------------------------


