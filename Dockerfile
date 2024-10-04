FROM mambaorg/micromamba

WORKDIR /home/$MAMBA_USER

# root level stuff
USER root
RUN apt-get update
RUN apt-get install -y git

RUN mkdir /home/data
RUN chown $MAMBA_USER /home/data

USER $MAMBA_USER

# get source code
RUN git clone https://github.com/medema-group/BiG-SCAPE.git -b feature/repo-updates

WORKDIR /home/$MAMBA_USER/BiG-SCAPE

# set up environment
RUN micromamba install \
    # channel
    -c conda-forge \
    # accept all prompts
    -y \
    # name of env
    -n base \
    # get python
    -f environment.yml && \
    # clean up cache
    micromamba clean --all --yes

RUN chmod 777 bigscape.py

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "./bigscape.py" ]
CMD [ "--help" ]
