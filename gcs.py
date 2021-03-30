import os
from google.cloud import storage

bucket_name = 'ntd-disease-simulator-data'
client = storage.Client()
bucket = client.bucket( bucket_name )

# gcs_path = relative from bucket name, so e.g. diseases/trachoma/data/...
def get_blob( gcs_path ):
    blob = bucket.blob( gcs_path )
    bytes = blob.download_as_bytes()
    return bytes

def blob_exists( gcs_path ):
    blob = bucket.blob( gcs_path )
    found = blob.exists()
    return found

def write_string_to_file( string, gcs_path ):
    blob = bucket.blob( gcs_path )
    blob.upload_from_string( string )

def download_blob_to_file( gcs_path, file_path ):
    blob = bucket.blob( gcs_path )
    blob.download_to_filename( file_path )

def upload_file_to_blob( file_path, gcs_path ):
    blob = bucket.blob( gcs_path )
    blob.upload_from_filename( file_path )
