FROM pshannon/trenaserver:latest
USER root
WORKDIR /home/trena
RUN mkdir data
RUN mkdir pkgs
COPY server.R /home/trena/work/
COPY data/* /home/trena/data/
COPY pkgs/* /home/trena/pkgs/
USER root
RUN mkdir sharedData
RUN chmod 777 sharedData
RUN chown -R trena:users /home/trena/work
RUN chown -R trena:users /home/trena/data
RUN chown -R trena:users /home/trena/sharedData
USER trena
WORKDIR /home/trena/work
RUN R CMD INSTALL /home/trena/pkgs/BrowserViz_2.0.1.tar.gz
RUN R CMD INSTALL /home/trena/pkgs/trena_1.1.12.tar.gz
RUN R CMD INSTALL /home/trena/pkgs/trenaViz_0.99.39.tar.gz

CMD ["R", "-f", "server.R"]
