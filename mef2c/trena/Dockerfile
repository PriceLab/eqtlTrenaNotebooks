FROM pshannon/trenaserver:latest
USER root
WORKDIR /home/trena
RUN mkdir data
RUN mkdir pkgs
RUN mkdir tmp

COPY server.R /home/trena/work/
COPY newServer.R /home/trena/work/
COPY data/* /home/trena/data/
COPY pkgs/* /home/trena/pkgs/
RUN cat myBashCustomizations >> /home/trena/.bashrc

USER root
RUN mkdir sharedData
RUN chmod 777 sharedData
RUN chown -R trena:users /home/trena/work
RUN chown -R trena:users /home/trena/pkgs
RUN chown -R trena:users /home/trena/data
RUN chown -R trena:users /home/trena/sharedData
RUN chown -R trena:users /home/trena/tmp

USER trena
WORKDIR /home/trena/work

RUN R CMD INSTALL  --no-test-load /home/trena/pkgs/BrowserViz_2.0.1.tar.gz
RUN R CMD INSTALL  --no-test-load /home/trena/pkgs/trena_1.1.12.tar.gz
RUN R CMD INSTALL  --no-test-load /home/trena/pkgs/trenaViz_0.99.39.tar.gz
RUN R CMD INSTALL  --no-test-load /home/trena/pkgs/SingleGeneData_0.99.2.tar.gz
RUN R CMD INSTALL  --no-test-load /home/trena/pkgs/MEF2C.data_0.99.3.tar.gz
RUN R CMD INSTALL  --no-test-load /home/trena/pkgs/SingleGeneAnalyzer_0.99.8.tar.gz
RUN R CMD INSTALL  --no-test-load /home/trena/pkgs/MotifDb_1.21.2.tar.gz

CMD ["R", "-f", "newServer.R"]
