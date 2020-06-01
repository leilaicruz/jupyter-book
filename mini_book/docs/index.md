# To reproduce an example of  Jupyter book, and you are a Windows user:

**Do this**: 


The following workflow should succeed using a miniconda powershell terminal on Windows 10:

`Open Anaconda Power Shell Promt`

0. `conda install -c eoas_ubc runjb`
1. `conda install git`
2. `git clone https://github.com/eoas-ubc/quantecon-mini-example.git`
3. `cd quantecon-mini-example`
4. `git checkout windows`
5. `conda env create -f environment_win.yml`
6. `conda activate wintest`
7. `cd mini_book`
8. `runjb docs`

After the build, view the html with (local deployment):

`start docs\_build\html\index.html`


