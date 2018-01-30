FROM pshannon/ipytrenaviz:v0.9.12
USER root
WORKDIR /home/jovyan
COPY myBashCustomizations /home/jovyan/
RUN cat myBashCustomizations >> /home/jovyan/.bashrc
#RUN mkdir /home/jovyan/.jupyter
COPY custom.css /home/jovyan/.jupyter/custom/

RUN mkdir /home/jovyan/mef2c
RUN mkdir /home/jovyan/mef2c/notebooks

RUN mkdir /home/jovyan/mef2c/hub
COPY hub/* /home/jovyan/mef2c/hub/
COPY notebooks/* /home/jovyan/mef2c/notebooks/

# ADD . /home/jovyan/mef2c/
RUN mkdir /home/jovyan/mef2c/privateData
RUN mkdir /home/jovyan/mef2c/sharedData

RUN chown -R jovyan:users /home/jovyan/.jupyter
RUN chown -R jovyan:users /home/jovyan/work
RUN chown -R jovyan:users /home/jovyan/mef2c

USER jovyan
WORKDIR /home/jovyan/mef2c
