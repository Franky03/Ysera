from flask import Flask, url_for, send_from_directory, request, render_template
import os
import requests
import json
from werkzeug.utils import secure_filename
from ysera import myfunction, ysera
from util import base64ToFile, fileToBase64
from flask_cors import CORS
import threading
import time

app: Flask = Flask(__name__)
CORS(app)

PROJECT_HOME = os.path.dirname(os.path.realpath(__file__))

def toYsera(params, pdbBase64, id):
    name = 'file_'+str(id)+'.pdb'
    base64ToFile(pdbBase64,PROJECT_HOME+"/temp/"+name)
    data = ysera(name, params)
    returnResult(id, data)


@app.route('/ysera', methods = ['POST'])
def getFile():
    
    content = request.get_json(silent=True)
    pdbBase64 = content["pdbBase64"]
    params = content["params"]
    id = content["id"]
    threading.Thread(target=toYsera, args=(params, pdbBase64, id)).start()

    return json.dumps("ok")


def returnResult(id, data):
    url = 'http://localhost:3030/calculate/'+str(id) + "/result"

    body = {
        "complete": "data:text/plain;base64, " + fileToBase64(PROJECT_HOME+"/output/file_"+str(id)+".pdb.txt"),
        "data" : data,
    }
    
    x = requests.post(url, data = body)
    print(x.text)

@app.route('/record', methods=['POST'])
def download():
    data = request.get_json()
    print(data) 
    
    return send_from_directory(directory='output', filename=data['filename'])
    

app.run(host='0.0.0.0', debug=True)
