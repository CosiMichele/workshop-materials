FROM stack0verflow/jupyterlab-pytorch

USER root

# Install RStudio
RUN apt update && apt install --yes gdebi-core && \
    wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2023.06.0-421-amd64.deb && \
    gdebi -n rstudio-server-2023.06.0-421-amd64.deb && \
    rm rstudio-server-2023.06.0-421-amd64.deb && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Shiny Server
RUN wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.20.1002-amd64.deb && \
    gdebi -n shiny-server-1.5.20.1002-amd64.deb && \
    rm shiny-server-1.5.20.1002-amd64.deb && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

USER jovyan
WORKDIR /home/jovyan

# Install Jupyter Lab Proxy extensions (cards in Launcher)
RUN mamba install --quiet --yes jupyter-rsession-proxy jupyter-server-proxy && \
    pip install jupyter-shiny-proxy
    
# Install and configure jupyter lab. 
COPY jupyter_notebook_config.json /opt/conda/etc/jupyter/jupyter_notebook_config.json

# Rebuild the Jupyter Lab with new tools
RUN jupyter lab build

EXPOSE 8888

COPY entry.sh /bin
COPY *ipynb /home/jovyan/
RUN mkdir -p /home/jovyan/.irods

ENTRYPOINT ["bash", "/bin/entry.sh"]