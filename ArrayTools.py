# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
#
# Ceci est un module du code d'Analyse des resultats des simulations qui contient un 
# ensemble de methods agissent sur les arrays.
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
import numpy as np

# ----------------------------------------------------------------------------------------------

# This method is used for extracting and reshaping a 2D array,
# x1-coordinate and y-coordinate values from a flattened table
# inputs:
#   extractValue: (float/int) value to be extracted
#   columIndexes: (integer32)(4) column indexes:
#		    0: index of the value to extract
#		    1: index of the x axis value
#		    2: index of the y axis value
#		    3: index of the value to plot
#   tableSize:    (set(integer32))(2) valueTable dimensions
#		    (nrows,ncolumns)
#   values:	    (float64)(N,M) value table to slice 
# outputs:
#   xAxisValue: (float64)(nrows) abscissa values
#   yAxisValue: (float64)(ncolumns) ordinate values
#   valueTable: (float64)(nrows,ncolumns) size of the table to be plotted
def extractArrayFromValueReshape(extractValue,columIndexes,tableSize,values):

  # extract test column
  testColumn = values[:,columIndexes[0]]
  # extract x values
  xAxisValue = np.unique(values[testColumn==extractValue,columIndexes[1]])
  # extract y values
  yAxisValue = np.unique(values[testColumn==extractValue,columIndexes[2]])
  # extract and reshape table value
  valueTable = np.transpose(\
  np.reshape(values[testColumn==extractValue,columIndexes[3]],tableSize)) 

  # return all values
  return valueTable,xAxisValue,yAxisValue

# ----------------------------------------------------------------------------------------------
