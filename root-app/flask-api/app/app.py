from flask import Flask, request, jsonify
from flask_cors import CORS #Flask-Cors
import pandas as pd
import copy

app = Flask(__name__)

CORS(app)

# Import functions 
from .engine import ( 
    get_unique_idCLs, 
    get_antibodies,
    plot_antibodies,
    get_celltypes,
    plot_celltypes,
    get_experiments,
    which_antibodies,
    which_celltypes,
    get_tissues,
)

@app.route('/idcls', methods=['GET'])
def idcls():
    uniqueIDCLs = get_unique_idCLs()
    result = {'idCLs': uniqueIDCLs}
    return result

@app.route('/findabs', methods=['POST'])
def findabs():
    res_dict = request.json
    background_idCLs = None
    idBTO_filter = None
    idExp_filter = None

    res_dict_keys = list(res_dict.keys())

    idCL_family_dict = res_dict['idCL']

    if 'background' in res_dict_keys:
        background_idCLs = res_dict['background']

    if 'idBTO' in res_dict_keys: 
        idBTO_filter = res_dict['idBTO']
    
    if 'idExp' in res_dict_keys:
        idExp_filter = res_dict['idExp']

    result_df = get_antibodies(node_fam_dict=idCL_family_dict,
                               custom_background_fam_dict=background_idCLs, 
                               idBTO=idBTO_filter, 
                               idExperiment=idExp_filter)
    
    if isinstance(result_df, str):
        return result_df
    else:
        result_json = result_df.to_json()
        return result_json

@app.route('/plotabs', methods=['POST'])
def plotabs():
    abs_idcls = request.json

    abs = abs_idcls['abs']
    idcls = abs_idcls['idcls']

    result_df = plot_antibodies(ab_ids=abs, 
                                id_CLs=idcls)
    result_json = result_df.to_json()

    return result_json

def nested_dicts(d):
    # https://stackoverflow.com/questions/54358643/python-nested-dictionaries-of-dataframes-convert-to-json
    for k, v in d.items():
        if isinstance(v, pd.DataFrame):
            d[k] = v.to_json()
        else:
            d[k] = nested_dicts(v)
    return d

@app.route('/findcelltypes', methods=['POST'])
def findcelltypes():
    res_dict = request.json
    idBTO_filter = None
    idExp_filter = None

    res_dict_keys = list(res_dict.keys())

    antibodies = res_dict['ab']

    if 'idBTO' in res_dict_keys: 
        idBTO_filter = res_dict['idBTO']
    
    if 'idExp' in res_dict_keys:
        idExp_filter = res_dict['idExp']

    result_dict = get_celltypes(ab_ids=antibodies, 
                                idBTO=idBTO_filter, 
                                idExperiment=idExp_filter)
    
    if isinstance(result_dict, str):
        return result_dict
    else:
        dict_of_dfs = nested_dicts(copy.deepcopy(result_dict))
        return dict_of_dfs

@app.route('/plotcelltypes', methods=['POST'])
def plotcelltypes():
    ab_idcls = request.json

    ab = ab_idcls['ab']
    idcls = ab_idcls['idcls']

    result_df = plot_celltypes(ab_id=ab, 
                               id_CLs=idcls)
    result_json = result_df.to_json()

    return result_json

@app.route('/findexperiments', methods=['POST'])
def findexperiments():
    res_dict = request.json
    idCL_filter = None
    idBTO_filter = None
    
    res_dict_keys = list(res_dict.keys())
    antibodies = res_dict['ab']

    if 'idCL' in res_dict_keys:
        idCL_filter = res_dict['idCL']

    if 'idBTO' in res_dict_keys: 
        idBTO_filter = res_dict['idBTO']

    result_df = get_experiments(ab_id=antibodies,
                                idCL=idCL_filter,
                                idBTO=idBTO_filter)
    
    if isinstance(result_df, str):
        return result_df
    else:
        result_json = result_df.to_json()
        return result_json

@app.route('/whichantibodies', methods=['POST'])
def whichantibodies():
    res_dict = request.json
    search_query = res_dict['search_query']

    result_df = which_antibodies(search_query=search_query)

    # If no antibodies from that search query were found in the database
    if isinstance(result_df, str):
        return result_df
    else:
        result_json = result_df.to_json()
        return result_json

@app.route('/whichcelltypes', methods=['POST'])
def whichcelltypes():
    res_dict = request.json
    search_query = res_dict['search_query']

    result_df = which_celltypes(search_query=search_query)

    # If no idCLs from that search query were found in the database
    if isinstance(result_df, str):
        return result_df
    else:
        result_json = result_df.to_json()
        return result_json

@app.route('/tissues', methods=['GET'])
def findtissues():
    tissues = get_tissues()
    tissues_json = tissues.to_json()
    return tissues_json

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000)