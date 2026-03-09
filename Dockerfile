# Build shiny app container for use with shinyproxy
FROM openanalytics/r-ver:4.5.2

# Install dependencies
RUN apt-get update && apt-get install -y zlib1g-dev

# Install packages for running shiny application

RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org/')"

RUN R -e "install.packages('shiny', repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('shinydashboard', repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('shinyjs', repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('ggplot2', repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('scales', repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('reshape', repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('bslib', repos='https://cloud.r-project.org/')"

RUN R -e "remotes::install_github(repo = 'SMAC-Group/gmwm', ref = 'gmwm2')"
RUN R -e "remotes::install_github(repo = 'SMAC-Group/wv')"

# copy the app to the image
RUN mkdir /root/gui4gmwm2
COPY R /root/gui4gmwm2/R

CMD ["R", "-e", "shiny::runApp('/root/gui4gmwm2/R', port = 3838, host='0.0.0.0')"]
