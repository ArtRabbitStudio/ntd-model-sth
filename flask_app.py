import json
import pandas
import hashlib

from gcs import blob_exists
from datetime import datetime
from flask import Flask, Response, request
from flask_cors import CORS, cross_origin
from sth_simulation.helsim_RUN import STH_Simulation

bucket_name = 'ntd-disease-simulator-data'
source_data_path_root = f"diseases/sth-roundworm/source-data"
source_data_gcs_path_root = f"/{bucket_name}/{source_data_path_root}"
output_data_path_root = f"diseases/sth-roundworm/data"
output_data_gcs_path_root = f"/{bucket_name}/{output_data_path_root}"
gs_prefix = "gs:/"
https_prefix = "https://storage.googleapis.com"

app = Flask(__name__)
cors = CORS( app, resources = { r"/run": { "origins": "*" } } ) # TODO FIXME to right origin
app.config[ 'CORS_HEADERS' ] = 'content-type'

@app.route('/')
def root():
    return Response( '👋\n', mimetype = 'text/plain' )

@app.route( '/run', methods = [ 'POST', 'OPTIONS' ] )
@cross_origin( origin = "*", headers = [ 'content-type' ] )
def run():

    # read in configuration from POST
    request_data_str = str( request.data, 'UTF-8' )
    request_hash = hashlib.sha256( request_data_str.encode( 'UTF-8' ) ).hexdigest()

    iu = request.json[ 'iu' ]
    region = iu[ 0:3 ]
    column_names = request.json[ 'mdaData' ][ 0 ]
    mda_data = request.json[ 'mdaData' ][ 1: ]

    GcsMDAFilePath = f"{gs_prefix}{output_data_gcs_path_root}/{region}/{request_hash}/InputMDA-{request_hash}.csv"
    HttpsMDAFilePath = f"{https_prefix}{output_data_gcs_path_root}/{region}/{request_hash}/InputMDA-{request_hash}.csv"

    try:

        pandas.DataFrame(
            request.json[ 'mdaData' ][ 1: ],
            columns = request.json[ 'mdaData' ][ 0 ]
        ).to_csv(
            GcsMDAFilePath,
            index = None
        )

    except Exception as e:

        return json.dumps( {
            'status': False,
            'msg': str( e )
        } )

    # ultimately going to come from input configuration
    GcsMDAFilePath = f"{gs_prefix}{source_data_gcs_path_root}/Asc_22609/STH_MDA_22609.csv"

    # going to come from cloud storage
    RkFilePath = f"{gs_prefix}{source_data_gcs_path_root}/Asc_22609/Input_Rk_Asc_22609.csv"
    InSimFilePath = f"{source_data_path_root}/Asc_22609/Asc_22609.p"

    # going to go to cloud storage
    PrevFileName = f"Output_Asc_22609.csv"

    # TODO FIXME sort this out
    PrevFilePath = f"{gs_prefix}{output_data_gcs_path_root}/{PrevFileName}"
    PrevBlobPath = f"{output_data_path_root}/{PrevFileName}"
    HttpsPrevFilePath = f"{https_prefix}{output_data_gcs_path_root}/{PrevFileName}"

    isNewSimulation = False

    if not blob_exists( PrevBlobPath ):

        isNewSimulation = True

        STH_Simulation(
            # comes from inside the python module
            paramFileName = "AscarisParameters_moderate.txt",
            demogName = "WHOGeneric",
            MDAFilePath = GcsMDAFilePath,
            PrevFilePath = PrevFilePath,
            RkFilePath = RkFilePath,
            nYears = 10,
            outputFrequency = 12,
            numReps = None,
            SaveOutput = False,
            OutSimFilePath = None,
            InSimFilePath = InSimFilePath,
            useCloudStorage = True
        )

    result = json.dumps( {
        'status': True,
        'scenarioHash': request_hash,
        'isNewSimulation': isNewSimulation
    } )

    try:
        return Response( result, mimetype = 'application/json; charset=UTF-8' )

    except Exception as e:
        return str( e )

if __name__ == '__main__':
    app.run( debug = False, host = '0.0.0.0' )
