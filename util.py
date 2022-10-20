import base64
import os

PROJECT_HOME = os.path.dirname(os.path.realpath(__file__))

def fileToBase64(path):
  with open(path, "rb") as f:
    data = f.read()
    
    decoded = base64.b64encode(data)
    return decoded.decode('utf-8')
  
def base64ToFile(string, path):
  w = open(path, 'w')
  
  splt = string.split(",")
  if (len(splt) > 1):
    string = splt[1]
  
  file = base64.b64decode(string)  
  w.write(file.decode('utf-8'))
  w.close()
  return path
  
# b64 = pdbToBase64(PROJECT_HOME + '/temp/teste.pdb')

# w = open(PROJECT_HOME + '/temp/base64.txt', 'w' )
# w.write(b64)
# w.close()

# rb64 = open(PROJECT_HOME + '/temp/base64.txt', 'r')
# data = rb64.read()

# base64ToPdb(data, PROJECT_HOME + '/temp/n.pdb')