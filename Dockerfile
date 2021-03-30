# vi: ft=Dockerfile
FROM python:3.8-buster

# get all the python libraries in early
RUN pip3 install numpy
RUN pip3 install pandas
RUN pip3 install joblib
RUN pip3 install flask
RUN pip3 install flask_cors
RUN pip3 install google-cloud-storage
RUN pip3 install fsspec
RUN pip3 install gcsfs
RUN pip3 install gunicorn

# put the app together
RUN mkdir -p /app
WORKDIR /app
COPY sth_simulation sth_simulation
COPY Pipfile Pipfile
COPY flask_app.py flask_app.py
COPY gcs.py gcs.py
COPY wsgi.py wsgi.py
COPY setup.py setup.py
RUN ls -l
RUN grep -vw sys < flask_app.py > removed-sys-flask_app.py && mv -f removed-sys-flask_app.py flask_app.py

# install the STH code
RUN pip3 install .

# port
EXPOSE 5000

# run the app in gunicorn via WSGI on port 5000
CMD gunicorn --workers 1 --timeout 600 --bind 0.0.0.0:5000 wsgi
