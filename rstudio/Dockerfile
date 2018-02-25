FROM rocker/rstudio:latest
USER root
RUN apt-get update
RUN apt-get install libxml2-dev -y
RUN apt-get -y install libmariadb-client-lgpl-dev -y
RUN apt-get install libpq-dev -y
RUN sudo apt-get install libpng-dev -y
RUN sudo apt-get install zlib1g-dev -y
RUN mkdir pkgs
COPY . /home/rstudio
COPY pkgs/* /home/rstudio/pkgs/
WORKDIR /home/rstudio
RUN mkdir /home/rstudio/library
RUN chmod 777 /home/rstudio/library
#RUN chown -R rstudio:users /home/

user rstudio
RUN cat /home/rstudio/myBashCustomizations >> /home/rstudio/.bashrc

RUN R -f installTrenaPrequisites.R

#RUN R CMD INSTALL  --no-test-load /home/rstudio/pkgs/BrowserViz_2.0.1.tar.gz
#RUN R CMD INSTALL  --no-test-load /home/rstudio/pkgs/trena_1.1.12.tar.gz
#RUN R CMD INSTALL  --no-test-load /home/rstudio/pkgs/trenaViz_0.99.39.tar.gz
#RUN R CMD INSTALL  --no-test-load /home/rstudio/pkgs/SingleGeneData_0.99.4.tar.gz
#RUN R CMD INSTALL  --no-test-load /home/rstudio/pkgs/MEF2C.data_0.99.15.tar.gz
#RUN R CMD INSTALL  --no-test-load /home/rstudio/pkgs/SingleGeneAnalyzer_0.99.18.tar.gz
#RUN R CMD INSTALL  --no-test-load /home/rstudio/pkgs/MotifDb_1.21.2.tar.gz
