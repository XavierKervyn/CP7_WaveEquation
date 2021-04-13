# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
#
# Ceci est un module du code d'Analyse des resultats des simulations qui contient un 
# ensemble de methods mathematiques.
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

# Cette fonction retourne a function pointer avec la fonction pour faire des
# derivees 1D
# inputs:
# derivativeType: (str) type de derivee 1D
# outputs:
# derivative1D: (funct) fonction pour calculer un derivee 1D
def select1DDerivative(derivativeType):
  # selectioner la derivee
  if(derivativeType=='forwardDerivative'):
    # retourner la derivee 1 order en avant
    return forwardDerivative1D1Order
  else:
    # warning
    print('Warning: 1D derivative type not found')
    # retourner
    return

# ----------------------------------------------------------------------------------------------

# Calculer les derivee en avant 1 ordre avec padding a 0
# inputs:
#   mesh:  (array(double))(NMesh) vecteur du maillage
#   value: (array(double))(NMesh) vecteur des valeurs
# outputs:
#  derivative: (array(double))(NMesh) vecteur des derivees
def forwardDerivative1D1Order(mesh,value):

  # initialiser array des derivee
  derivative = np.zeros((len(mesh)),dtype=np.float64)

  # boucle sur les valeurs de la mesh
  for nodeId,node in enumerate(mesh[1:]):
    # calculer la derivee
    derivative[nodeId+1] = (value[nodeId+1]-value[nodeId])/(node-mesh[nodeId])
  # retourner les valeurs
  return derivative

# ----------------------------------------------------------------------------------------------

# Cette fonction retourne a function pointer avec la fonction pour faire des
# interpolations
# inputs:
# interpolationType: (str) type d'interpolation
# outputs:
# derivative1D: (funct) fonction pour calculer une interpolation
def selectInterpolation(interpolationType):
  # selectioner la derivee
  if(derivativeType=='linear1D'):
    # retourner l'interpolation lineaire 1D
    return linearInterpolation1D
  else:
    # warning
    print('Warning: interpolation type not found')
    # retourner
    return

# ----------------------------------------------------------------------------------------------

# Claculer un interpolation lineaire 1D entre deux points
# inputs
#   position:  (float64) position d'interpolation
#   nodes:     (float64)(2) 0: position 1iere node
#			      1: position 2ieme node
#   nodeValue: (float64)(2) 0: valeur 1iere node
#			      1: valeur 2ieme node
# outputs
#   value: (float64) interpolated value 
def linearInterpolation1D(position,nodes,nodeVaues):

  # calculer l'interpolation
  return (nodeValue[0]+((nodeValue[1]-nodeValue[0])*(position-nodes[0]))/(nodes[1]-nodes[0]))

# ----------------------------------------------------------------------------------------------

# Calculer le maillage 1D pour interpolations 1D lineaires uniformes
# inputs:
#   normalisedPosition: (float64) position normalise des points
#			d'interpolations dans la mesh
#   meshNodes:		(array(float64))(NPoints) array containing
#			the original mesh node positions  
# outputs:
#   interpolationMesh:  (array(float64))(NPoints-1) array containing
#   			uniform interpolation nodes
def computeUniformMeshLinearInterpolation1D(normalisedPosition,meshNodes):

  # initialise the interpolation mesh
  interpolationMesh = np.zeros((meshNodes.size[0]-1,1),dtype=meshNodes.dtype)

  # loop on the original mesh nodes
  for nodeId,node in enumerate(meshNodes[0:-1]):
    # compute the interpolation mesh
    interpolationMesh[nodeId] = node+normalisedPosition*(node[nodeId+1]-node)

  # return the iterpolation mesh
  return interpolationMesh

# ----------------------------------------------------------------------------------------------

# Retourner la valeur du premier pique en utilisant la method find_peaks
def first_peak_value(data):

  # importer le module signal du scipy
  from scipy import signal as signal

  # obtenir les piques et leur indexes  
  peak_values,peak_ids = peak_values_indexes(data)
  # retourner la veleur du premier pique
  return peak_values[0]

# ----------------------------------------------------------------------------------------------

# Retourner les valeurs et les indexes des piques en utilisant la method find_peaks
def peak_values_indexes(data):

  # importer le module signal du scipy
  from scipy import signal as signal

  # extraire les indexes des piques
  peak_ids = signal.find_peaks(data)
  # extraire les valuers des piques
  peak_values = data[peak_ids[0]]
  # retourner les valeurs des piques et leur indexes
  return peak_values,peak_ids[0]

# ----------------------------------------------------------------------------------------------

