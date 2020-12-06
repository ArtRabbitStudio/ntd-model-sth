import json
import sys
from google.cloud import storage

bucket_name = 'ntd-disease-simulator-data'
client = storage.Client()
bucket = client.bucket( bucket_name )

iu_ids = json.load( open( sys.argv[ 1 ] ) )

for iu_id in iu_ids:
    country = iu_id[:3]
    padded_id = iu_id[3:].zfill(5)
    path = f"diseases/{sys.argv[2]}/source-data/{country}/{iu_id}/{sys.argv[3]}_{iu_id}.p"
    blob = bucket.blob( path )
    if not blob.exists():
        print( f"==> {iu_id} not found", flush=True )
    else:
        print( f"    {iu_id} found OK", flush=True )

