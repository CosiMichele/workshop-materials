FROM cosimichele/jupyter-amber:gpu-250625

USER jovyan
WORKDIR /home/jovyan

RUN pip install git+https://github.com/CharlieLaughton/Alphafix.git
RUN pip install git+https://github.com/CompBioAsia/CBAtools.git
RUN pip install --upgrade parmed

# Rebuild the Jupyter Lab with new tools
RUN jupyter lab build

EXPOSE 8888

#COPY entry.sh /bin

#ENTRYPOINT ["bash", "/bin/entry.sh"]
