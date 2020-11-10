import json
import pandas
import hashlib

from gcs import blob_exists, write_string_to_file
from datetime import datetime
from flask import Flask, Response, request
from flask_cors import CORS, cross_origin
from sth_simulation.helsim_RUN import STH_Simulation

bucket_name = 'ntd-disease-simulator-data'
gs_prefix = "gs:/"
https_prefix = "https://storage.googleapis.com"

parameter_file_names = {
    'sth-roundworm': "AscarisParameters_moderate.txt"
}

file_name_disease_abbreviations = {
    'sth-roundworm': "Asc"
}

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
    request_hash = hashlib.sha256( request_data_str.encode( 'UTF-8' ) ).hexdigest()[ 0:24 ]


    # whip out necessary vars
    iu = request.json[ 'iu' ]
    country = iu[ 0:3 ]
    iu_id = iu[ 3: ]
    column_names = request.json[ 'mdaData' ][ 0 ]
    mda_data = request.json[ 'mdaData' ][ 1: ]
    numReps = 200 if request.json[ 'runs' ] > 200 else request.json[ 'runs' ]

    disease = request.json[ 'disease' ]
    paramFileName = parameter_file_names[ disease ]
    file_abbrev = file_name_disease_abbreviations[ disease ]


    # set up all the file paths
    source_data_path_root = f"diseases/{disease}/source-data"
    source_data_gcs_path_root = f"/{bucket_name}/{source_data_path_root}"

    output_data_path_root = f"diseases/{disease}/data"
    output_data_gcs_path_root = f"/{bucket_name}/{output_data_path_root}"

    OutputDirectoryPath = f"{output_data_path_root}/{country}/{iu}/{request_hash}"
    OutputDirectoryGsPath = f"/{bucket_name}/{OutputDirectoryPath}"

    MDAFilePath = f"{OutputDirectoryPath}/InputMDA-{request_hash}.csv"
    GcsMDAFilePath = f"{gs_prefix}/{bucket_name}/{MDAFilePath}"

    GcsRkFilePath = f"{gs_prefix}{source_data_gcs_path_root}/{country}/{iu}/Input_Rk_{file_abbrev}_{iu}.csv"
    HttpsRkFilePath = f"{https_prefix}{source_data_gcs_path_root}/{country}/{iu}/Input_Rk_{file_abbrev}_{iu}.csv"
    InSimFilePath = f"{source_data_path_root}/{country}/{iu}/{file_abbrev}_{iu}.p"

    PrevFileName = f"OutputPrevKKSAC-{file_abbrev}-{iu}-{request_hash}.csv"

    PrevBlobPath = f"{OutputDirectoryPath}/{PrevFileName}"
    GcsPrevFilePath = f"{gs_prefix}{OutputDirectoryGsPath}/{PrevFileName}"
    HttpsPrevFilePath = f"{https_prefix}{OutputDirectoryGsPath}/{PrevFileName}"

    PrevSummaryFileName = PrevFileName[:-4] + "-summary.json"
    GcsPrevSummaryFilePath = f"{OutputDirectoryPath}/{PrevSummaryFileName}"
    HttpsPrevSummaryFilePath = f"{https_prefix}{OutputDirectoryGsPath}/{PrevSummaryFileName}"

    HistoricalPrevFileName = f"PrevKKSAC{file_abbrev}_{iu}.csv"
    GcsHistoricalPrevFilePath = f"{gs_prefix}{source_data_gcs_path_root}/{country}/{iu}/{HistoricalPrevFileName}"
    HttpsHistoricalPrevFilePath = f"{https_prefix}{source_data_gcs_path_root}/{country}/{iu}/{HistoricalPrevFileName}"

    HistoricalPrevSummaryFileName = f"HistoricalPrev-{iu}-{request_hash}-summary.json"
    GcsHistoricalPrevSummaryFilePath = f"{OutputDirectoryPath}/{HistoricalPrevSummaryFileName}"
    HttpsHistoricalPrevSummaryFilePath = f"{https_prefix}{OutputDirectoryGsPath}/{HistoricalPrevSummaryFileName}"

    Result = {
        'status': True,
        'isNewSimulation': False,
        'futureDataUrl': HttpsPrevFilePath,
        'futureSummaryUrl': HttpsPrevSummaryFilePath,
        'historicalDataUrl': HttpsHistoricalPrevFilePath,
        'historicalSummaryUrl': HttpsHistoricalPrevSummaryFilePath
    }


    # convert the incoming scenario mdaData to a CSV and write it to GCS
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


    # run the scenario, if its output hasn't already been written to cloud storage
    if not blob_exists( PrevBlobPath ):

        # we're about to kick off a new simulation
        Result[ 'isNewSimulation' ] = True

        STH_Simulation(
            paramFileName = paramFileName, # comes from inside the python module
            demogName = "WHOGeneric", # standard for STH
            MDAFilePath = GcsMDAFilePath,
            PrevFilePath = GcsPrevFilePath,
            RkFilePath = GcsRkFilePath,
            nYears = 12,
            outputFrequency = 6, # restrict the number of columns in the CSV output
            numReps = numReps,
            SaveOutput = False,
            OutSimFilePath = None,
            InSimFilePath = InSimFilePath,
            useCloudStorage = True
        )

        # summarize generated future prevalence data (predictions)
        future_prevalence = pandas.read_csv( GcsPrevFilePath )
        future_summary = pandas.DataFrame( {
            'median': future_prevalence.iloc[:, 2:].median(),
            'lower': future_prevalence.iloc[:, 2:].quantile(0.05),
            'upper': future_prevalence.iloc[:, 2:].quantile(0.95)
        }).to_json()
        write_string_to_file( future_summary, GcsPrevSummaryFilePath )

        # summarize historica prevalence data
        historical_prevalence = pandas.read_csv( GcsHistoricalPrevFilePath )
        historical_summary = pandas.DataFrame( {
            'median': historical_prevalence.iloc[:, 2:].median(),
            'lower': historical_prevalence.iloc[:, 2:].quantile(0.05),
            'upper': historical_prevalence.iloc[:, 2:].quantile(0.95)
        }).to_json()
        write_string_to_file( historical_summary, GcsHistoricalPrevSummaryFilePath )


    try:

        # snag the output for sending to browser now
        output_result_json = json.dumps( Result )

        # save result to file for JS to hit next time
        ResultJsonFilePath = f"{OutputDirectoryPath}/{file_abbrev}-{iu}-{request_hash}-info.json"
        Result[ 'isNewSimulation' ] = False # because reading from static file means it's not new
        write_string_to_file( json.dumps( Result ), ResultJsonFilePath )

        return Response( output_result_json, mimetype = 'application/json; charset=UTF-8' )

    except Exception as e:

        return json.dumps( {
            'status': False,
            'msg': str( e )
        } )


if __name__ == '__main__':
    app.run( debug = False, host = '0.0.0.0' )
