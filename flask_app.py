import json

from gcs import blob_exists
from datetime import datetime
from flask import Flask, Response
from sth_simulation.helsim_RUN import STH_Simulation

bucket_name = 'ntd-disease-simulator-data'
source_data_path_root = f"diseases/sth/source-data"
source_data_gcs_path_root = f"/{bucket_name}/{source_data_path_root}"
output_data_path_root = f"diseases/sth/data"
output_data_gcs_path_root = f"/{bucket_name}/{output_data_path_root}"
gs_prefix = "gs:/"
https_prefix = "https://storage.googleapis.com"

app = Flask(__name__)

@app.route('/')
def hello():

    # ultimately going to come from input configuration
    MDAFilePath = f"{gs_prefix}{source_data_gcs_path_root}/Asc_22609/STH_MDA_22609.csv"


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
            demogName = "KenyaKDHS",
            MDAFilePath = MDAFilePath,
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
        'PrevFilePath': HttpsPrevFilePath,
        'IsNewSimulation': isNewSimulation
    } )

    try:
        return Response( result, mimetype = 'application/json; charset=UTF-8' )

    except Exception as e:
        return str( e )

if __name__ == '__main__':
    app.run( debug = False, host = '0.0.0.0' )
