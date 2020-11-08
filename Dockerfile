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

# put the app together
RUN mkdir -p /app
WORKDIR /app
COPY sth_simulation sth_simulation
COPY Pipfile Pipfile
COPY flask_app.py flask_app.py
COPY gcs.py gcs.py
COPY setup.py setup.py
COPY files/Asc_22609/ files/Asc_22609/
RUN ls -l

# install the STH code
RUN pip3 install .

# port
EXPOSE 5000

ENTRYPOINT [ "/usr/local/bin/python" ]
CMD [ "flask_app.py" ]
