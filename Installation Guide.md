# BiG-SCAPE installation

Although each library could be installed on its own, the use of a virtual
environment is highly recommended. Here is a quick guide of BiG-SCAPE 
installation using Miniconda

* Install [Miniconda](https://conda.io/miniconda.html). This will install Python 2 as default for all new conda environments. You'll need to re-log for the changes to go into effect.

```
> wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
> chmod u+x Miniconda2-latest-Linux-x86_64.sh
> ./Miniconda2-latest-Linux-x86_64.sh
```

* Create new environmnent:

```
> conda create --name bigscape
```

* Activate new environmnent

```
> source activate bigscape
```

* Install packages:

```
> conda install numpy scipy
> conda install -c bioconda hmmer biopython mafft
```

* The Affinity Propagation algorithm used currently is [pySAPC](https://pypi.python.org/pypi/pysapc/1.1.0). If you are using a Mac: `conda install -c https://conda.anaconda.org/bioinfocao pysapc`, otherwise: `pip install pysapc`

(pysapc also installs pandas, cython, pytz, six and python-dateutil)
