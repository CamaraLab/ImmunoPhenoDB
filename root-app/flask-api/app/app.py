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
    which_experiments,
    get_tissues,
    downsample_reference_table,
    database_statistics,
    antibody_panel_reference_table,
    get_antibodies_web,
    plot_antibodies_web,
    get_celltypes_web,
    plot_celltypes_web,
    get_experiments
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
        elif isinstance(v, dict):  # If the value is a nested dictionary, recursively apply the function
            d[k] = nested_dicts(v)
        elif isinstance(v, str):  # If the value is a string, leave it as is
            d[k] = v
        else:
            raise TypeError(f"Unsupported type for key {k}: {type(v)}")
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
    
    print("findcelltypes result_dict:", result_dict, type(result_dict))
    
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
    
@app.route("/whichexperiments", methods=['POST'])
def whichexperiments():
    res_dict = request.json
    search_query = res_dict['search_query']

    result_df = which_experiments(search_query=search_query)
    # If no experiments from that search query
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

@app.route('/stveareference', methods=['POST'])
def stveareference():
    res_dict = request.json
    print("stveareferernce dict:\n", res_dict)
    antibody_pairs = res_dict["antibody_pairs"]
    idBTO = res_dict["idBTO"]
    idExperiment = res_dict["idExperiment"]
    parse_option = res_dict["parse_option"]
    # pairwise_threshold = res_dict["pairwise_threshold"]
    # na_threshold = res_dict["na_threshold"]
    population_size = res_dict["population_size"]

    result_df = downsample_reference_table(antibody_pairs=antibody_pairs, 
                                           idBTO=idBTO,
                                           idExperiment=idExperiment,
                                           parse_option=parse_option,
                                           population_size=population_size)
    result_json = result_df.to_json()
    return result_json

@app.route('/databasestatistics', methods=['GET'])
def databasestatistics():
    stats_json = database_statistics()
    return stats_json

@app.route('/antibodypanelreference', methods=['POST'])
def antibodypanelreference():
    res_dict = request.json
    print("antibodypanelreference dict:", res_dict)
    target_idcl = res_dict["target_idcl"]
    target_idbto = res_dict["target_idbto"]
    background_idcl = res_dict["background_idcl"]
    background_idbto = res_dict["background_idbto"]
    experiment = res_dict["experiment"]
    seed = res_dict["seed"]

    result_df = antibody_panel_reference_table(unique_target_family_idCLs=target_idcl,
                                               final_target_idBTOs=target_idbto,
                                               modified_background_family_idCLs=background_idcl,
                                               final_background_idBTOs=background_idbto,
                                               experiment=experiment,
                                               seed=seed)
    result_json = result_df.to_json()
    return result_json

# Website functions
@app.route('/findabsweb', methods=['POST'])
def findabsweb():
    # JSON payload needs to contain idCLs and idBTO
    res_dict = request.json

    idCL_query_list = res_dict['idCL']
    idBTO_query_list = res_dict['idBTO']

    result_df = get_antibodies_web(idCLs=idCL_query_list,
                                   idBTO=idBTO_query_list)
    
    # Check if an error string was produced
    if isinstance(result_df, str):
        return result_df
    else:
        # If no error, convert the dataframe into a JSON to send over network
        result_json = result_df.to_json()
        return result_json

@app.route('/plotabsweb', methods=['POST'])
def plotabsweb():
    # JSON payload needs ab_ids and id_CLs. ab_ids is ideally retrieved from findabsweb
    # from a separate API call in the frontend
    res_dict = request.json

    ab_ids = res_dict['abs']
    id_CLs = res_dict['idcls']
    
    result_df = plot_antibodies_web(ab_ids=ab_ids,
                                    id_CLs=id_CLs)
    
    return result_df.to_json()

@app.route('/findcelltypesweb', methods=['POST'])
def findcelltypesweb():
    # JSON payload needs ab_ids and idBTO
    res_dict = request.json

    ab_query_list = res_dict['ab']
    idBTO_query_list = res_dict['idBTO']

    result_dict = get_celltypes_web(ab_ids=ab_query_list,
                                  idBTO=idBTO_query_list)
    
    print("findcelltypesweb result_dict:", result_dict, type(result_dict))
    
    # Check if an error string was produced
    if isinstance(result_dict, str):
        return result_dict
    else:
        dict_of_dfs = nested_dicts(copy.deepcopy(result_dict))
        return dict_of_dfs


@app.route('/plotcelltypesweb', methods=['POST'])
def plotcelltypesweb():
    # JSON payload needs one ab_id and id_CLs
    res_dict = request.json

    ab_id = res_dict['ab']
    id_CLs = res_dict['idcls']
    
    result_df = plot_celltypes_web(ab_id=ab_id,
                                   id_CLs=id_CLs)
    
    return result_df.to_json()

@app.route('/experiments', methods=['GET'])
def experiments():
    experiments_df = get_experiments()
    return experiments_df.to_json()

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000)