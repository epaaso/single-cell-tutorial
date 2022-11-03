sudo docker run --interactive --tty --name sc_anal --publish 8888-8892:8888-8892 --volume $HOME/2021-SC-HCA-LATAM/CONTAINER:/root/host_home --workdir /root/host_home/ leanderd/single-cell-analysis:210114 /bin/bash

jupyter trust ~/host_home/single-cell-tutorial/latest_notebook/Case-study_Mouse-intestinal-epithelium_2101.ipynb

pip install -U "ipython>=7.20"
