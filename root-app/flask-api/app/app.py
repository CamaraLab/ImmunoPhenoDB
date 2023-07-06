from flask import Flask, request, jsonify
from flask_cors import CORS
import pandas as pd

app = Flask(__name__)

CORS(app)

# Import functions 
from .engine import get_unique_idCLs, get_Abs

@app.route('/idCLs', methods=['GET'])
def idCLs():
    uniqueIDCLs = get_unique_idCLs()
    result = {'idCLs': uniqueIDCLs}
    print(result)
    return result

@app.route('/abLMM', methods=['POST'])
def runLMM():
    idCL_family_dict = request.json
    result_df = get_Abs(idCL_family_dict)
    result_json = result_df.to_json()
    print(result_json)
    return result_json

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000)
    # app.run(host="localhost", port=5000)