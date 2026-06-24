FROM cosimichele/jupyter-rfdiff:gpu-260615

USER jovyan
WORKDIR /home/jovyan

RUN pip install "jax[cuda12_pip]==0.4.34"

# Rebuild the Jupyter Lab with new tools
RUN jupyter lab build

EXPOSE 8888

