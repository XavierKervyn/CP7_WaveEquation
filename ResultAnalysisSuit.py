# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
#
# Ceci est un module du code d'Analyse des resultats des simulations qui a comme objectifs
# de calculer certain grandeurs physiques depuis les resultats des simulations.
# Le code d'analyse est base sur le standard python3 (on recommende python3.6 ou superieur).
# Pour installer les librairies on recommende d'utiliser pip3. Si pip3 n'est pas installe,
# il est possible de suivre la procedure d'un des links suivants:
#   https://linuxconfig.org/how-to-install-pip-on-ubuntu-18-04-bionic-beaver
#   https://linuxize.com/post/how-to-install-pip-on-ubuntu-18.04/
# 
# Ensuite, il faut installer les librairies: 
#   numpy 
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

# creer la class ResultAnalysisSuite
class ResultAnalysisSuit():

  # Definir le constructeur de la class PlotResultats
  # inputs:
  #   fileNameList:	       (list(string)) list des nomes du fichier avec les valeurs
  #   values:		       (list(array))(Nfiles,Nx,Ny) tableau des valeus
  #   parameters:	       (dict) dictionaire contenant les parametres pour l'analyse
  #				      des resultats
  #   figureParameters:	       (dict) parametres pour faire des figures
  #   valuesDelimiter:	       (string) signe pour delimiter les colonnes
  #					dans le fichier: default = ' '
  # outputs:
  # constructeur de la class
  def __init__(self,fileNameList=[],values=np.array([]),parameters={},\
  figureParameters={},valuesDelimiter=' '):

    # initialiser la liste des valeurs
    self.__values = []
    # verifier s'il faut lire ou sauvegarder les resultats
    if(len(fileNameList)!=0):
      # boucle sur les nomes des fichiers a lire
      for fileName in fileNameList:  
        # charger le tableau des resultats dans un numpy array
        self.__values.append(np.loadtxt(fileName,dtype=np.float64,delimiter=valuesDelimiter))
    else:
      # sauvegarder les donnees
      self.__values = values

    # sauvegarder les parametres pour analyser les resultats, parametres:
    # distance: (dict) dictionary containing parameters for computing phase space distances
    # 		indexes:       (list(list(int)))(NDistance,NIndex) list of list
    #					containing the value indexes to be used for 
    #					computing each distance
    #		normalisation: (list(list(double)))(NDistance,NIndex) list of list
    #					containing the normalisation to be used for 
    #					computing each distance for each index
    #		meshIndex:     (int) index of the mesh for producing plots
    #	
    self.parameters = parameters

    # sauvegarder les parametres pour generer des figures, parametres:
    #   plotData:       list(dict) contient tous informations pour faire des plots
    #                          chaque cles est une figure differents. Chaque cle
    #                                contient un dictionaire avec les donnees:
    #                                NRows:   (int) numero de lignes de subplots 
    #				     NCols:   (int) numero de colonnes de subplots
    #				     NPlots:  (list(int)) numero de lignes dans le meme
    #						    subplot pour chaque subplot
    #				     indexes: (list(list(list(int))))(NSubplots,NPlots,2)
    #						    liste d'indexes des colonnes 
    #						    a imprimer. La premier liste definit
    #						    la deuxieme les lignes sur le meme subplot
    #						    la troisieme les index sous la forme
    #						    [xId,yId] where xId is the index 
    #						    of des abscisses et yId celui des ordonnees
    #				     title:   (list(string))(NSubplots) liste de titre pour
    #						    chauque subplot
    #				     xLables: (list(string))(NSubplots) liste de string containant
    #						    nome des x-axes pour chaque subplot
    #				     yLables: (list(string))(NSubplots) liste de string containant
    #						    nome des y-axes pour chaque subplot
    #				     normalisedPlot: (list(string))(NSubplots) liste de bool
    #						    si vrais le plot est normalise
    #				     axis: (list(string))(NSubplots) axis aspect ratio
    #				     legend: (list(list(string)))(NSubplots,NPlots)
    # 				     		    add a legend to subplot
    #				     legendLocation: (list(string))(NSubplots,NPlots)
    # 				     		     location of the legend
    #
    #   fontSize:            (integer) figure font size
    #   lineStyle:           (list(string)) list of line stile for plotting
    #					    one for each value array
    #   lineWidth:           (integer) width of the line to plot
    #   marker:              (list(string)) list of types of marker to use:
    #					    one for each value array 
    #   markerSize:          (integer) size of the marker to use for plotting
    #   plotType:            (str) type of plot to be done, available:
    #			           1) linear, 2) loglog
    #   AnalyticalSolutions: (function)(NSubplots) fonction for solution comparison
    self.figureParameters = figureParameters

    # cette dictionaire contient les resultats de l'analyse
    self.__results = {} 

# ----------------------------------------------------------------------------------------------

    # destructeur de la class
    def __del__(self):

      # imprimer le nome de la class
      print(self.__class__.__name__ ,'is deleted')

# ----------------------------------------------------------------------------------------------

  # cette methode rendre le tableau des valeurs
  # inputs:
  # outputs:
  #   values: (array)(Nx,Ny) result table
  def getValues(self):

    # rendre le tableau des valeurs
    return self.__values

# ----------------------------------------------------------------------------------------------

  # cette methode copie le tableau des valeurs
  # inputs:
  #   values: (array)(Nx,Ny) result table
  # outputs:
  def setValues(self,values):

    # copier le tableau des valeurs
    self.__values = np.copy(values)

# ----------------------------------------------------------------------------------------------

  # cette methode rendre un reultats en fonction de la clee
  # inputs:
  #   key:    (string) clee du resultat
  # outputs:
  #   result: (array)(Nx,Ny) result table
  def getValues(self,key):

    # check if key in results
    if(key in self.__results):
      # rendre le resultat
      return self.__results[key]
    else:
      # print error
      print('key ',key,' not found!')
      # return None
      return

# ----------------------------------------------------------------------------------------------

  # cette methode copie un reultats en fonction de la clee
  # inputs:
  #   key:    (string) clee du resultat
  #   result: (array) tableau du resultat
  # outputs:
  def setValues(self,key,result):

    # copier le tableau des valeurs
    self.__results[key] = np.copy(result)

# ----------------------------------------------------------------------------------------------

  # cette methode est un interface pour executer differents type d'analyse
  # inputs:
  # outputs:
  def doResultAnalysis(self):

    # verifier quelle type d'analyse est demandee
    if('distance' in self.parameters):
      # calculer la distance entre deux solutions
      self.computeDistance()
    # verifier quelle type d'analyse est demandee
    if('computeCompareDerivative1D' in self.parameters):
      # executer comparaison entre un derivee et un valeur
      self.computeAndCompare1DDerivatives()
    if('findMaximaAndVelocity2D' in self.parameters):
      # compute maxima and maxima velocity in 2D tables
      self.computeMaximaAndVelocity2D()

# ----------------------------------------------------------------------------------------------

  # calculer la distance entre deux solutions
  # tous exposants sont calcules par rapport
  # premier tableau des valeurs
  # inputs:
  # outputs:
  def computeDistance(self):

    # importer le module pour faire des figures
    import PlotResults

    # Verifier la presence de au moins 2 tableau des valeurs
    if(len(self.__values)<2):
      # envoyer un warning
      print('Not enough values for computing distance between solutions')
      # sortir par la fonction
      return

    # initialiser le dictionare des resultats
    self.__results['distance'] = []

    # calculer le nombre de distance
    NDistances = len(self.parameters['distance']['indexes'])
    # calculer le nombre de valeurs
    NValues = len(self.__values)-1
       
    # boucle sur les values de 1 a NValues-1
    for valueId,value in enumerate(self.__values[1:]):
      # initialiser tableau des distances NMesh x NDistances
      distance = np.zeros((value.shape[0],NDistances),dtype=np.float64)
      # boucle sur les distances a calculer
      for indexListId,indexList in enumerate(self.parameters['distance']['indexes']):
        # boucle sur la liste des indexes
        for indexId,index in enumerate(indexList):
          # calculer les caree des distances par rapport au tableau 0
          distance[:,indexListId] = distance[:,indexListId] + \
          (value[:,index]-self.__values[0][:,index])*\
          (value[:,index]-self.__values[0][:,index])*\
          self.parameters['distance']['normalisation'][indexListId][indexId]*\
          self.parameters['distance']['normalisation'][indexListId][indexId]
        # calculer la racine caree
        distance[:,indexListId] = np.sqrt(distance[:,indexListId])
        # sauvegarder les valeur dans results
        self.__results['distance'].append(distance)

    # creer le dictionaire des parametres
    plotParameters = self.figureParameters.copy()
    # remove all other plot data
    plotParameters['plotData']=[self.figureParameters['plotData'][0]]
    # reset the number of lines per subplot
    plotParameters['plotData'][0]['NPlots'] = list(NValues*\
    np.ones((len(self.__results['distance'][0])),dtype=np.int32))
    # reset the set of indexes
    plotParameters['plotData'][0]['indexes'] = []
    # initialiser le tableau pour faire des figures
    figureValue = np.zeros((self.__values[0].shape[0],\
    NDistances*NValues+1),dtype=np.float64)
    # copier la mesh
    figureValue[:,0] = self.__values[0][:,self.parameters['distance']['meshIndex']]
    # boucle sur les distances
    for distanceArrayId,distanceArray in enumerate(self.__results['distance']):
      # initialise a sequence of indexes per subplots
      subplotIndexes = []
      # boucle sur les distances
      for indexId in range(distanceArray.shape[1]):
        # copier les resultats dans le tableau
        figureValue[:,NValues*indexId+distanceArrayId+1] = distanceArray[:,indexId]
        # append the new index
        subplotIndexes.append([0,NValues*indexId+distanceArrayId+1])
        # append the subplot list of indexes
        plotParameters['plotData'][0]['indexes'].append(subplotIndexes)

    # do figures
    # initialiser un objet pour generer des figures
    plotResult = PlotResults.PlotResults(values=[figureValue],\
    figureParameters=plotParameters)
    # fare la figure
    plotResult.SimplePlotFigures()
    
# ----------------------------------------------------------------------------------------------


  # Cette methode calcule la dervivee 1D pour un mesh et un index et le fait un plot avec un valeur
  # de comparaison
  # Si le vectteur des derivees est plus petit que celui des valeurs, les indices qui ne sont
  # pas utilies sont mis a zero.
  # inputs:
  # outputs:
  def computeAndCompare1DDerivatives(self):

    # importer le module pour faire des figures
    import PlotResults
    # importer des methodes mathematique
    from MathematicalTools import select1DDerivative

    # selectioner la methode pour les derivee
    derivative1D = select1DDerivative(\
    self.parameters['computeCompareDerivative1D']['derivativeType'])
  
    # calculer le numero des valeurs
    NValues = len(self.__values)
    # initialiser le dictionare des resultats
    self.__results['derivatives1D'] = np.zeros((self.__values[0].shape[0],\
    NValues*3),dtype=np.float64)

    # initialiser la figure
    # creer le dictionaire des parametres
    plotParameters = self.figureParameters.copy()
    # remove all other plot data
    plotParameters['plotData']=[self.figureParameters['plotData'][0]]
    # reset the number of lines per subplot
    plotParameters['plotData'][0]['NPlots'] = [NValues*2]
    # reset the set of indexes
    plotParameters['plotData'][0]['indexes'] = []
  
    # initialiser un subplot unique
    subplotIndexes = []
    # boucle sur le numero des valeurs
    for valueId,value in enumerate(self.__values):
      # sauvegarder la mesh
      self.__results['derivatives1D'][:,valueId] = \
      value[:,self.parameters['computeCompareDerivative1D']['meshIndex']]
      # appliquer la derivee
      self.__results['derivatives1D'][:,valueId+1] = derivative1D(\
      value[:,self.parameters['computeCompareDerivative1D']['meshIndex']],\
      value[:,self.parameters['computeCompareDerivative1D']['index']])
      # sauvegarder la valeur de comparaison
      self.__results['derivatives1D'][:,valueId+2] = \
      value[:,self.parameters['computeCompareDerivative1D']['comparisonIndex']]
      # sauvegarder les indices pour le plot
      plotParameters['plotData'][0]['indexes'].append([[valueId,valueId+1],[valueId,valueId+2]])

    # do figures
    # initialiser un objet pour generer des figures
    plotResult = PlotResults.PlotResults(values=[self.__results['derivatives1D']],\
    figureParameters=plotParameters)
    # fare la figure
    plotResult.SimplePlotFigures()

# ----------------------------------------------------------------------------------------------

  # Cette methode est utilise pour calculer les maximum d'une solution 2D
  # et ces vitesses. La solution et les maillages doivent etre donnes sous
  # la forme d'un tableau 2D ayant comme premier ligne la premiere maillage
  # et comme premiere colonne le deuxieme maillage. 
  def computeMaximaAndVelocity2D(self):

    # importer le module pour faire des figures
    import scipy.signal as signal
    import PlotResults

    # extrare les maillages et le tableau
    value = self.__values[0]
    xAxis = value[0,1:] # extraire 1iemere maillage
    yAxis = value[1:,0] # extraire 2xieme maillage
    values = np.array(value[1:,1:]) # couper le tableau
    # transposer le maillage
    if(self.parameters['findMaximaAndVelocity2D']['axis']=='x'):
      values = np.transpose(values)       
    # initialiser le tableau des resultats
    self.__results['findMaxima2D'] = np.zeros((len(xAxis),3),dtype=np.float64)
    self.__results['findVelocity2D'] = np.zeros((len(xAxis)-1,2),dtype=np.float64)
    # sauvegarder maillaige pour plot
    self.__results['findMaxima2D'][:,0] = xAxis[:]
    self.__results['findVelocity2D'][:,0] = np.array([0.5*(xAxis[elementId]+element) \
    for elementId,element in enumerate(xAxis[1:])])

    # generatre kernel for smoothing
    if('kernelSize' in self.parameters['findMaximaAndVelocity2D']):
      smoothingKernel = np.ones((\
      self.parameters['findMaximaAndVelocity2D']['kernelSize']),dtype=np.float64)/ \
      float(self.parameters['findMaximaAndVelocity2D']['kernelSize'])

    # boucle sur le premier index
    for columnId,column in enumerate(values):
      localValues = np.zeros((len(column)),dtype=np.float64)
      localValues[:] = column[:]
      # check if only branch has to be used
      if(self.parameters['findMaximaAndVelocity2D']['onlyPositiveNegative']=='positive'):
        localValues[column<0] = 0.0 
      elif(self.parameters['findMaximaAndVelocity2D']['onlyPositiveNegative']=='negative'):
        localValues[column>0] = 0.0
      # chercher les maxima en utilisant find peaks
      #localMaximaId = signal.argrelextrema(localValues,np.greater)
      localMaximaId = signal.find_peaks(localValues)
      if(len(localMaximaId[0])>0):
        localMaximaId = localMaximaId[0][0]
        if(localMaximaId>0):
          # sauvgarder valeur du 1iere maximum
          self.__results['findMaxima2D'][columnId,1] = column[localMaximaId]
          self.__results['findMaxima2D'][columnId,2] = yAxis[localMaximaId]

    # if present smooth the velocity
    if('kernelSize' in self.parameters['findMaximaAndVelocity2D']):
      # smooth the velocities
      self.__results['findMaxima2D'][:,2] = np.convolve(\
      self.__results['findMaxima2D'][:,2],\
      smoothingKernel,mode='same')

      # calculer la vitesse
      for timeId,time in enumerate(self.__results['findMaxima2D'][1:,2]):
        self.__results['findVelocity2D'][timeId,1] = \
        (xAxis[timeId+1]-xAxis[timeId])/ \
        (time-self.__results['findMaxima2D'][timeId,2])

    # initialiser la figure
    # creer le dictionaire des parametres
    plotParameters0 = self.figureParameters.copy()
    plotParameters1 = self.figureParameters.copy()
    # remove all other plot data
    plotParameters0['plotData']=[self.figureParameters['plotData'][0]]
    plotParameters1['plotData']=[self.figureParameters['plotData'][1]]
    # set first plot
    # subplot
    plotParameters0['plotData'][0]['NRows'] = 1
    plotParameters0['plotData'][0]['NCols'] = 1
    # reset the number of lines per subplot
    plotParameters0['plotData'][0]['NPlots'] = [1]
    # reset the set of indexes
    plotParameters0['plotData'][0]['indexes'] = [[[0,1]]]
    # set second plot
    # subplot
    plotParameters1['plotData'][0]['NRows'] = 1
    plotParameters1['plotData'][0]['NCols'] = 1
    # reset the number of lines per subplot
    plotParameters1['plotData'][0]['NPlots'] = [1]
    # reset the set of indexes
    plotParameters1['plotData'][0]['indexes'] = [[[0,1]]]
      
    # do figures
    # initialiser un objet pour generer des figures
    plotResult = PlotResults.PlotResults(values=[self.__results['findMaxima2D']],\
    figureParameters=plotParameters0)
    # fare la figure
    plotResult.SimplePlotFigures()
    # initialiser un objet pour generer des figures
    plotResult = PlotResults.PlotResults(values=[self.__results['findVelocity2D']],\
    figureParameters=plotParameters1)
    # fare la figure
    plotResult.SimplePlotFigures()


# ----------------------------------------------------------------------------------------------

