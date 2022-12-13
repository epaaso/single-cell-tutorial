# Instrucciones mejoras del README.md

Primero debes de correr el contenedor con su volumen externo para que permanezan los cambios

```bash
sudo docker run --interactive --tty --name sc_anal --publish 8888-8892:8888-8892 --volume $HOME/2021-SC-HCA-LATAM/CONTAINER:/root/host_home --workdir /root/host_home/ leanderd/single-cell-analysis:211119 /bin/bash
```

Solo una vez debes de confiar en el notebook,
para poder correr todos los comandos y asgurate de que la version de ipython sea 7.20:

```bash
jupyter trust ~/host_home/single-cell-tutorial/latest_notebook/Case-study_Mouse-intestinal-epithelium_2101.ipynb
pip install -U "ipython>=7.20"
```

Los cambios al notebook ya estan integrados en commits propios. Para poder empezar a trabajr con el notebook
debes de correr `jl`.
