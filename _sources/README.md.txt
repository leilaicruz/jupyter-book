[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/leilaicruz/jupyter-book/gh-pages)


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

In order to publish the Book online it is necessary to add to the build.yml:

- A public and private key just dedicated to the Jupyter Book repo 
- The content of the private key should be copied to the `Secrets` tab in the Settings of the repo. 
    - Create an SSH key with sufficient access privileges. For security reasons, don't use your personal SSH key but set up a dedicated one for use in GitHub Actions. See below for a few hints if you are unsure about this step.
    - Make sure you don't have a passphrase set on the private key.
    - In your repository, go to the Settings > Secrets menu and create a new secret. In this example, we'll call it SSH_PRIVATE_KEY. Put the contents of the private SSH key file into the contents field.
    - This key should start with `-----BEGIN ... PRIVATE KEY-----`, consist of many lines and ends with `-----END ... PRIVATE KEY-----`.
- SSH private key format
    - If the private key is not in the PEM format, you will see an Error loading key `"(stdin)": invalid format message`.

    - Use `ssh-keygen -p -f path/to/your/key -m pem` to convert your key file to PEM, but be sure to make a backup of the file first 

My build.yml file looks like:

```yaml
name: Test Build
on:
  push:
    branches-ignore:
      - master
jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout 🛎️
        uses: actions/checkout@v2
        with:
          persist-credentials: false

      - name: Setup Miniconda
        uses: goanpeca/setup-miniconda@v1
        with:
          auto-update-conda: true
          auto-activate-base: false
          miniconda-version: 'latest'
          python-version: 3.7
          environment-file: environment_win.yml
          activate-environment: wintest

      - name: Install jupyter_book
        shell: bash -l {0}
        run:  pip install jupyter-book --pre 


      - name: Build Jupyter Book
        shell: bash -l {0}
        run: jb build  mini_book/docs/

      - name: Install SSH Client 🔑
        uses: webfactory/ssh-agent@v0.2.0
        with:
          ssh-private-key: ${{ secrets.DEPLOY_KEY }}
        


      - name: Deploy 🚀
        uses: JamesIves/github-pages-deploy-action@releases/v3
        with:
          SSH: true
          BRANCH: gh-pages
          FOLDER: mini_book/docs/_build/html/
          
```
- Add an empty file `.nojekyll` to the root of the repo to prevent github on deploying a jekyll website.

