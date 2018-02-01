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

* Create new environmnent (you can choose Python 3 using ``python=3`` at the end of the next command):

```
> conda create --name bigscape
```

* Activate new environmnent

```
> source activate bigscape
```

* Install packages:

```
> conda install numpy scipy scikit-learn
> conda install -c bioconda hmmer biopython mafft fasttree
> conda install -c anaconda networkx
```
