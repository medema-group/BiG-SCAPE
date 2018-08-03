FROM conda/miniconda3

RUN apt-get update && apt-get install -y git wget

RUN conda create --name bigscape
#RUN source activate bigscape
RUN [ "/bin/bash", "-c","source activate bigscape"]
RUN conda install -y \
 	numpy \  
	scipy \ 
	scikit-learn 

RUN conda install -c bioconda hmmer biopython fasttree
RUN conda install -c anaconda networkx

WORKDIR /usr/src
## Cloning BiG-SCAPE
RUN git clone https://git.wur.nl/medema-group/BiG-SCAPE.git

## geting Pfam
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
#ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz.
RUN gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm && mv Pfam-A.* /usr/src/BiG-SCAPE/.

RUN chmod +x /usr/src/BiG-SCAPE/*py  
RUN chmod 777 /home  
ENV PATH /usr/src/BiG-SCAPE:$PATH
USER 1000:1000
RUN mkdir /home/input /home/output
WORKDIR /home
ENTRYPOINT ["bigscape.py"]
CMD ["--help"]

