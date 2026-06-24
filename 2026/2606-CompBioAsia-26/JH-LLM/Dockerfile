FROM cosimichele/jupyter-llm:260612

RUN npm install -g @anthropic-ai/claude-code

RUN npm install -g @musistudio/claude-code-router

COPY start.sh /home/jovyan/start.sh 
RUN sudo chmod 777 /home/jovyan/start.sh

# Rebuild the Jupyter Lab with new tools
RUN jupyter lab build

EXPOSE 8888
