# HERMES Burst Search Tools

The HERMES Burst Search Tools is a software collection of tools for finding transients
in data from the HERMES nanosatellite constellation. This software was developed in a collaboration between
the Italian National Institute for Astrophysics and the Italian Space Agency.

It is composed of a main library, called `hbstools`, and a collection of scripts to interface with it.
HERMES Burst Search Tools uses a cutting edge changepoint detection algorithm called Poisson FOCuS for discovering new astrophysical events. See:

* _Ward, K., Dilillo, G., Eckley, I., & Fearnhead, P. (2023). Poisson-FOCuS: An efficient online method for detecting count bursts with application to gamma ray burst detection. Journal of the American Statistical Association, 1-13._
* _Dilillo, G., Ward, K., Eckley, I. A., Fearnhead, P., Crupi, R., Evangelista, Y., Vaccjo. A., & Fiore, F. (2023). Gamma-ray burst detection with Poisson-FOCuS and other trigger algorithms. arXiv preprint arXiv:2312.08817._

----

## Setup:

### Downloading 

To download this repository, click on the "Code" green button at top right, then click on "Download ZIP".
Unzip the archive where best suited. This will become your installation folder.

Otherwise, if you have `git` installed, `cd` to the installation folder and run `git clone xxx`.


----
### Installing with Anaconda

Supposing you downloaded hbstools to `/path_to/hbstools`, run:

1. `cd /path_to/hbstools`
2. `conda create -n hbstools python=3.11 poetry`
3`conda activate hbstools`
4`poetry install`

Done!

> ❓ **Are you using Anaconda?**: 
> Try calling `which python` in your terminal. If you are using anaconda, 
> it should return something like: `/Users/username/!!anaconda3!!/envs/hbstools/bin/python`.

To test that everything is working try launching the Python console and `import hbstools`, or
try `mercury --help`.

> ❗ **Remember to always activate your environment**, otherwise you won't be able to use hbstools or mercury.
----
### Installing with venv

If you are installing with a virtual environment, `cd` to the folder where you downloaded hbstools, then run:

1. `python3 -m venv .`
2. `source ./bin/activate` or just `./bin/activate` if you are on Windows
3. `pip install poetry`
4. `poetry install`

----
### Uninstalling
If you are not using [poetry](https://python-poetry.org/) for other projects, you can run the following command, else refer to [this link](https://python-poetry.org/docs/managing-environments/#:~:text=Deleting%20the%20environments,-Finally%2C%20you%20can&text=You%20can%20delete%20more%20than%20one%20environment%20at%20a%20time.&text=Use%20the%20%2D%2Dall%20option%20to%20delete%20all%20virtual%20environments%20at%20once.&text=If%20you%20remove%20the%20currently,it%20will%20be%20automatically%20deactivated.).

1. `poetry env remove --all`

If you installed with Anaconda, also run `conda env remove -n hbstools`.

----

# Mercury

![mercury](assets/mercury-gif/mercury.gif)

Mercury is the first interface to `hbstools`. 
It is a command line tool to search for gamma-ray bursts and other astronomical transients. 

> 💅 **Mercury is best rendered on modern terminal applications.**
> If you are working on windows, we suggest using the [new Windows terminal](https://apps.microsoft.com/store/detail/windows-terminal/9N0DX20HK701).

### Searching for GRBs
A basic usage involves getting into the folder containing the data you want to analyze and running:

```mercury search .```

By default, only the input directory and its subdirectories are searched for data, but you can search deeper using the
`--reclim` option of `mercury search`.

> ❗ **To get help with mercury run `mercury --help`.**
> To get help on a particular command, such as `search`, you call `mercury search --help`.

### Configuration files

Mercury requires a configuration to work. Even when you run  `mercury search .` a default configuration is loaded.
You can get a configuration stub using:

```mercury drop .```

This will create a commented configuration file which you can modify with a text editor. 
To run `mercury search` with configuration `myconfig.yml`, run `mercury search -c myconfig.yml .`.
You can get a less verbose configuration file using `mercury --quiet drop .`

### Results

By default results are saved in the input directory in FITS format.
You can change the output file using `mercury search . -o myresults-filename.fits`

----

## Demo dataset

We have uploaded a demo dataset [online](https://drive.google.com/file/d/1kC473-QQsLWrClxKRHT8JJCIJr_KO_4_/view?usp=sharing).
Download the archive and unzip it:

1. `conda activate hbstools`
2. `cd /path_to/demodataset`
3. `mercury search .`

This supposing you installed `hbstools` with Anaconda. If you didn't, go to your virtual environment's bin directory and activate it with `source activate`.
