# HERMES Burst Search Tools
The HERMES Burst Search Tools is a software collection of tools for finding transients
in data from the HERMES nanosatellite constellation. This software was developed in a collaboration between
the Italian National Institute for Astrophysics and the Italian Space Agency.

It is composed of a main library, called `hbstools`, and a collection of scripts to interface with it.
HERMES Burst Search Tools uses a changepoint detection algorithm called Poisson FOCuS to seach for astrophysical events, see:

* _Ward, K., Dilillo, G., Eckley, I., & Fearnhead, P. (2023). Poisson-FOCuS: An efficient online method for detecting count bursts with application to gamma ray burst detection. Journal of the American Statistical Association, 1-13._
* _Dilillo, G., Ward, K., Eckley, I. A., Fearnhead, P., Crupi, R., Evangelista, Y., Vacchi, A., & Fiore, F. (2023). Gamma-ray burst detection with Poisson-FOCuS and other trigger algorithms. arXiv preprint arXiv:2312.08817._


# Setup
### Download
First download this repository. 
If you have `git` installed, you can do so running from terminal:

`git clone https://github.com/peppedilillo/hbstools.git`

This will create a directory called `hbstools` with the content of this repository.
Otherwise, you can click on the "Code" green button at top right, then click on "Download ZIP" and unzip the archive where best suited.

### Installing with Anaconda
Supposing you downloaded this repository to `/path_to/hbstools`, run:

1. `cd /path_to/hbstools`
2. `conda create -n hbstools python=3.11 poetry`
3. `conda activate hbstools`
4. `poetry install`

Done!

To test that everything is working try launching a Python console and `import hbstools`, or
run `mercury --help` from the terminal.

> ❗ **Remember to activate your environment**, otherwise you won't be able to use hbstools or mercury.

### Installing with venv
If you are installing with a virtual environment, move with terminal to the installation folder then run:

1. `python3 -m venv .`
2. `source ./bin/activate` or just `./bin/activate` if you are on Windows
3. `pip install poetry`
4. `poetry install`

### Uninstalling
From terminal run:

1. `poetry env remove --all`
2. `conda env remove -n hbstools`

If you are using poetry for other projects refer to [this link](https://python-poetry.org/docs/managing-environments/#deleting-the-environments) instead.
The latter step is only required if you installed with Anaconda.


# Mercury
![mercury](assets/mercury-gif/mercury.gif)

Mercury is the first interface to `hbstools`. 
It is a command line tool to automatically search for gamma-ray bursts and other astronomical transients. 

> 💅 **Mercury is best rendered on modern terminal applications.**
> If you are working on windows, we suggest using the [new Windows terminal](https://apps.microsoft.com/store/detail/windows-terminal/9N0DX20HK701).

### Searching for GRBs
A basic usage involves getting into the folder containing the data you want to analyze and running:

```mercury search .```

By default, only the input directory and its subdirectories are searched for data, but you can search deeper using the
`--reclim` option of `mercury search`.

### Configuration files
Mercury requires a configuration to work. Even when you run  `mercury search .` a default configuration is loaded.
You can get a configuration stub using:

```mercury drop .```

This will create a commented configuration file which you can later modify with a text editor. 
To run with a custom configuration use the `-c` option flag, e.g. `mercury search -c myconfig.yml`.

### Results
By default, results are saved in the input directory in FITS format.
Using `mercury search . -o myresults-filename.fits` you will change the output file to `myresults-filename.fits`.

> ❗ **To get help with mercury run `mercury --help`.**
> To get help on a particular command, such as `search`, you call `mercury search --help`.

### Demo dataset
We have uploaded a demo dataset [online](https://drive.google.com/file/d/1kC473-QQsLWrClxKRHT8JJCIJr_KO_4_/view?usp=sharing).
Download the archive and unzip it then, from terminal:

1. `conda activate hbstools`
2. `cd /path_to/demodataset`
3. `mercury search .`

This supposing you installed `hbstools` with Anaconda. 
If you didn't, activate your local environment instead.
