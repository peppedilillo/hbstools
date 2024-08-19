from glob import glob
import hashlib
from math import log10
from pathlib import Path
from typing import Sequence
from uuid import uuid4
import warnings

from astropy.io import fits
import click
import pandas as pd
from rich.console import Console
from schema import And  # type: ignore[import-untyped]
from schema import Optional  # type: ignore[import-untyped]
from schema import Schema  # type: ignore[import-untyped]
from schema import SchemaError  # type: ignore[import-untyped]
from schema import Use  # type: ignore[import-untyped]
import yaml
from yaml import dump as write_yaml
from yaml import safe_load as read_yaml
from yaml import YAMLError

import hbstools as hbs
from hbstools.data import map_event_to_files
from hbstools.types import Dataset
from hbstools.types import Event

LOGO = """
                        
  .uef^"                
:d88E                   
`888E            uL     
 888E .z8k   .ue888Nc.. 
 888E~?888L d88E`"888E` 
 888E  888E 888E  888E  
 888E  888E 888E  888E  
 888E  888E 888E  888E  
 888E  888E 888& .888E  
m888N= 888> *888" 888&  
 `Y"   888   `"   "888E 
      J88"  .dWi   `88E 
      @%    4888~  J8%  
    :"       ^"===*"`   

"""


def init_console(with_logo) -> Console:
    """initializes console, optionally with a logo"""

    def print_logo(c: Console):
        """Prints a logo header."""
        c.print("[white]" + LOGO, highlight=False)
        return c

    console = Console(width=80, highlight=False)
    if with_logo:
        print_logo(console)
    return console


DEFAULT_CONFIG = """
# The `binning` keyword sets the light curve bin step duration.
# It is expressed in units of seconds. 
# Must be greater than zero.
binning: 0.1

# Data are filtered over one customizable energy band.
# The `energy_lims` keyword sets the low and high energy limit for this band.
# Energy limits are expressed in units of keV.
# The first value (the low energy limit) must be non negative and smaller 
# than the second value (the high energy limit).
energy_lims:
- 20.0 # low energy limit
- 300.0 # high energy limit

# The `skip` keyword determines for how long the algorithm doesn't work
# after a trigger is encountered. This prevents re-trigger from long bursts
# with complex light curves or high variability.
# It is expressed in units of bin-steps. This means that if `binning` is set
# to 0.1 and `skip` is set to 3000, the algorithm won't work for the 5 minutes
# which follow a trigger. Passed the skip time, the algorithm is restarted.
# Must be a non-negative integer.
skip: 3000


algorithm_params:
 # Different algorithms support different parameters. 
 # The algorithm we are going to use depends on the parameters you'll select, 
 # according to this table:
 # 
 #  +---------------+---------------+------------+----------+-------+
 #  |  parameters   | Python PF+DES | Python BFT | C PF+SES | C BFT |
 #  +---------------+---------------+------------+----------+-------+
 #  | threshold_std |       v       |     v      |    v     |   v   |
 #  |    mu_min     |       v       |     v      |    v     |   v   |
 #  |     alpha     |       v       |     v      |    v     |   v   |
 #  |     beta      |       v       |     v      |          |       |
 #  |       m       |       v       |     v      |    v     |   v   |
 #  |     t_max     |       v       |     v      |          |       |
 #  |     sleep     |       v       |     v      |    v     |   v   |
 #  |   majority    |               |     v      |          |   v   |
 #  +---------------+---------------+------------+----------+-------+
 #
 # NOTE: a successful execution is not guaranteed if your parameters selection
 # does not meet one of the columns in the above table.
 # NOTE: `t_max` and `beta` are commented by default.

 # The `threshold` keyword sets the confidence level for detection.
 # The higher the threshold, the less false positive but the more false negative.
 # The threshold is expressed in units of standard deviations.
 # Must be greater than 0.
 threshold_std: 4.5

 # The `mu_min` key set a focus parameter for killing old changepoints, which 
 # most likely will never result in a trigger. 
 # Not used if equal 1.0, keep it below 1.5.
 # Must not be smaller than 1.
 mu_min: 1.1
 
 # The `alpha` keyword sets the background smoothing characteristic time tau,
 # > tau = (binning/0.005).
 # Must be greater than 0.
 alpha: 0.005

 # The `beta` keyword sets the trend component for background estimate.
 # The algorithm may become unstable when `beta` is set: keep it off (0.0) 
 # unless you have good reasons to turn it on.
 # Must be non-negative.
 # beta: 0.0

 # The `m` keyword will prevent the most recent observed count to be used for 
 # background estimation. This prevents background estimate to be "polluted"
 # by real transients. 
 # It is expressed in units of bin-steps. This means that if `binning` is set
 # to 0.1 and `m` is set to 40, the algorithm won't work use the latest 4.0 s
 # of data in background estimate. 
 # Must be a positive integer.
 m: 40

 # The `t_max` parameter tells the algorithm to kill old changepoint.
 # It is a good idea to keep it equal to `m`. Also expressed as a bin-step.
 # Must be a positive integer.
 # t_max: 40

 # The algorithms stays idle for a while before startin its operations.
 # This help forming a good estimate of the background. 
 # The `sleep` parameters sets for how long this idle period lasts.
 # In particular, testing for anomalies will start after `m + sleep` bin-steps.
 # The `sleep` parameter is expressed as a bin-step.
 # Must be a non-negative integer.
 sleep: 120
 
 # Running one of the BFT algorithms you may select how many detectors must
 # be simultaneously over threshold for a trigger to pass through.
 # This is done specifying the `majority` parameter.
 # Must be an integer comprised between 1 and 4 (included).
 majority: 3
 """


config_schema = Schema(
    {
        "binning": And(
            Use(float),
            lambda b: b > 0,
            error="`binning` must be greater than zero.",
        ),
        "energy_lims": And(
            Use(tuple[float, float]),
            lambda lim: 0 <= lim[0] < lim[1],
        ),
        "skip": And(
            lambda s: isinstance(s, int) & (s >= 0),
            error="`skip` must be a non-negative integer",
        ),
        "algorithm_params": {
            "threshold_std": And(
                Use(float),
                lambda t: t > 0,
                error="`threshold` must be positive.",
            ),
            "alpha": And(
                Use(float),
                lambda a: a > 0,
                error="`alpha` must be positive",
            ),
            Optional("beta"): And(
                Use(float),
                lambda b: b >= 0,
                error="`beta` must be non-negative",
            ),
            "mu_min": And(
                Use(float),
                lambda m: m >= 1,
                error="`mu_min` must be equal or greater than one",
            ),
            "m": And(
                lambda m: isinstance(m, int) & (not m < 1),
                error="`m` must be a positive integer.",
            ),
            "sleep": And(
                lambda sleep: isinstance(sleep, int) & (not sleep < 0),
                error="`sleep` must be a non-negative integer.",
            ),
            Optional("t_max"): And(
                lambda t_max: isinstance(t_max, int) & (not t_max < 1),
                error="`t_max` must be a positive integer.",
            ),
            "majority": And(
                lambda maj: isinstance(maj, int) & (not maj < 1) & (not maj > 4),
                error="`majority` must be an integer between 1 and 4 (included).",
            ),
        },
    },
    ignore_extra_keys=True,
)


def crawler(
    directory: Path | str,
    patterns: Sequence[str],
    recursion_limit: int = 1,
) -> set[tuple[Path]]:
    """Goes through directories looking for subdirs containing a specific set of files,
    each matching unix-style patterns only once.
    If `recursion_limit=0`  we only check the present folder, ignoring its subdirectories.
    Returns a set (unordered) of tuples. Each tuple is composed of `n` path-objects,
    where `n` equals the number of patterns. Each path-object points to a matching file in
    a subdirectory. Partial matches are discarded. An error is raised if more than one
    file matches a pattern.
    """

    def get_subdirs(d: Path) -> list[Path]:
        return [f for f in d.iterdir() if f.is_dir()]

    def directory_contains(d: Path, ps) -> tuple:
        matches = {p: [d / f for f in glob(p, root_dir=d)] for p in ps}
        # if we do not get a full match we return without errors
        if any([not matches[p] for p in ps]):
            return tuple()
        # raise an error if we got an ambiguous match
        if not all(len(matches[p]) == 1 for p in ps):
            raise click.FileError(
                f"Directory {d} contains multiple files pattern matching against "
                f"`{', '.join([p for p in ps if len(matches[p]) != 1])}`."
            )
        return tuple(m for p in ps for m in matches[p])

    def check_subdirs(
        subdirs: list[Path], ps, rlim: int, recursion_acc: int, acc: set[tuple]
    ):
        if not subdirs:
            return []
        car, *cdr = subdirs
        check_dir(car, ps, rlim, recursion_acc + 1, acc)
        check_subdirs(cdr, ps, rlim, recursion_acc, acc)

    def check_dir(d: Path, ps, rlim, recursion_acc: int, acc: set[tuple]) -> set[tuple]:
        if recursion_acc == rlim + 1:
            return acc
        if matches := directory_contains(d, ps):
            acc.add(matches)
        check_subdirs(get_subdirs(d), ps, rlim, recursion_acc, acc)
        return acc

    return check_dir(Path(directory), patterns, recursion_limit, 0, set())


def unused_path(file: Path, num: int = 1, isdir: bool = False) -> Path:
    """Return an unused path with same stem prefix as `file` and incremental suffix.
    if `num` is set to 1, returns `path.fits`, `path-1.fits`, `path-2.fits` etc.
    if `num` is set to 0, return `path.fits`, `path-0.fits`, `path-1.fits`..`"""
    if (not file.is_file() and not isdir) or (not file.is_dir() and isdir):
        return file
    parts = file.stem.split("-")
    *tail, head = parts
    stem = "-".join(tail) if head.isdigit() else "-".join(parts)
    fname = f"{stem}-{num}{file.suffix}"
    return unused_path(file.parent / fname, num=num + 1, isdir=isdir)


def fmt_filename(filename: str | Path) -> str:
    """Rich formatting for filenames."""
    return f"'[b]{filename}[/]'"


@click.group()
@click.option(
    "--quiet",
    "-q",
    type=click.BOOL,
    is_flag=True,
    help="Write less words.",
)
@click.version_option(package_name="hbstools")
@click.help_option()
@click.pass_context
def cli(ctx: click.Context, quiet: bool):
    """This is mercury a command line interface to hbstools for searching
    transients in data from the HERMES spacecrafts."""
    ctx.obj = {
        "quiet": quiet,
    }
    return


def search_validate_config(
    ctx: click.Context, param: click.Option, config_path: Path | None
) -> dict:
    """Validates user configuration option, and make sure it matches one of
    the available trigger algorithms."""

    if config_path is None:
        config = read_yaml(DEFAULT_CONFIG)
    else:
        with open(config_path, "r") as stream:
            config = read_yaml(stream)

    try:
        validated_config = config_schema.validate(config)
    except YAMLError:
        raise click.BadParameter("Cannot parse YAML configuration file.")
    except SchemaError as error:
        raise click.BadParameter(f"Wrong inputs in config file: \n - {error}")

    # check if we have an algorithm compatible with the configuration.
    algorithm_params = validated_config["algorithm_params"]
    try:
        algorithm_class = hbs.trigger.trigger_match(algorithm_params)
    except ValueError:
        raise click.BadParameter("Cannot find an algorithm matching the configuration.")

    # stores the algorithm's name and the configuration path to display them later.
    ctx.obj["search_algoname"] = str(algorithm_class(**algorithm_params))
    ctx.obj["search_config"] = None if config_path is None else config_path
    return validated_config


DEFAULT_CATALOG_NAME = "mercury-results.fits"


def write_catalog(events: list[Event], filepath: Path | str):
    """Save an event list to a single fits file"""
    content = pd.DataFrame(events).to_records(index=False)
    fits.writeto(filename=filepath, data=content, overwrite=True)


DEFAULT_LIBRARY_NAME = "mercury-results/"
INDEX_FILENAME = ".mercury-index.yaml"


def write_library(events: list[Event], dataset: Dataset, dir_path: Path):
    """
    Writes events to multiple fits files and store an index of them which
    can be used to merge the results to the dataset.

    :param events:
    :param dataset:
    :param dir_path: shall point to a directory. both the events and the index
    are given default names.
    """

    def write_src(event: Event, filepath: Path, header: fits.Header | None = None):
        """A helper for writing an event's source output file"""
        content = pd.DataFrame([event])[["start", "end"]].to_records(index=False)
        fits.writeto(filename=filepath, data=content, header=header)

    def write_bkg(e: Event, filepath: Path, header: fits.Header | None = None):
        """A helper for writing an event's background output file"""
        content = pd.DataFrame(
            {
                "bkg_start": [e.bkg_pre_start, e.bkg_post_start],
                "bkg_end": [e.bkg_post_start, e.bkg_post_end],
            }
        ).to_records(index=False)
        fits.writeto(filename=filepath, data=content, header=header)

    index = {"uuid": uuid4().hex, "mappings": (fmap := {})}
    pad = int(log10(len(events))) + 1  # for filename padding
    for n, (event, (_, gti_path)) in enumerate(
        map_event_to_files(events, dataset).items()
    ):
        header = hbs.io.read_gti_header(gti_path)
        write_src(event, src_path := dir_path / f"event-src-{n:0{pad}}.fits", header)
        write_bkg(event, bkg_path := dir_path / f"event-bkg-{n:0{pad}}.fits", header)
        fmap[src_path.name] = {"root": str(gti_path.parent.absolute()), "type": "src"}
        fmap[bkg_path.name] = {"root": str(gti_path.parent.absolute()), "type": "bkg"}

    with open(dir_path / INDEX_FILENAME, "w") as f:
        yaml.dump(index, f)


@cli.command()
@click.argument(
    "input_dirs",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    nargs=-1,
)
@click.option(
    "-o",
    "--output",
    "output",
    type=click.Path(dir_okay=True, path_type=Path),
    default=Path("."),
    help="Path to a file or a directory where results are saved. "
    "If not provided, the output are saved in the input directory. ",
)
@click.option(
    "-c",
    "--config",
    "configuration",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    callback=search_validate_config,
    help="Path to a YAML configuration file. "
    "If not provided, values from a default configuration are used.",
)
@click.option(
    "-m",
    "--mode",
    "mode",
    type=click.Choice(["catalog", "library"]),
    default="catalog",
    help="Formats the output. In `catalog` mode a single fits file is created."
    "In `library` mode, a folder or fits file is generated",
)
@click.option(
    "--reclim",
    "recursion_limit",
    type=click.INT,
    default=1,
    help="Mercury uses a filesystem crawler to locate data in subdirectories."
    "This parameter sets a limit to recursion depth, how deep subfolders are searched.",
)
@click.pass_context
def search(
    ctx: click.Context,
    input_dirs: tuple[Path],
    output: Path,
    configuration: dict,
    mode: str,
    recursion_limit: int,
):
    """Searches transients from data in the input directories.
    The algorithm used for this search is called Poisson-FOCuS (Ward, 2022; Dilillo, 2024).
    The search is configurable using a YAML configuration, see mercury's 'drop' command.
    """
    # fmt: off
    console = init_console(with_logo=True)
    console.print(f"Welcome, this is [bold]mercury.search[/].\n")
    config_file = "default" if (f := ctx.obj["search_config"]) is None else f
    console.log(f"Loaded {fmt_filename(config_file)} configuration.")
    console.log(f"Algorithm {ctx.obj['search_algoname']} matches the configuration.")

    patterns = ["*_lv1_cl.evt", "*gti.fits"]
    data_paths = {
        subdir for directory in input_dirs
        for subdir in crawler(directory, patterns, recursion_limit)
    }
    if not data_paths:
        console.print("\nFound no data. Exiting.\n")
        return
    console.log(f"Found {len(data_paths)} data folder{'' if len(data_paths) == 1 else 's'}.")

    dataset = hbs.data.catalog(data_paths)
    events = hbs.search(dataset, configuration, console=console)
    if not events:
        console.print("\nNo results to save. Exiting.\n")
        return

    if mode == "catalog":
        filepath = unused_path(output if not output.is_dir() else output / DEFAULT_CATALOG_NAME)
        console.log(f"Writing to {fmt_filename(filepath)}.")
        write_catalog(events, filepath)
    elif mode == "library":
        dirpath = unused_path(output / DEFAULT_LIBRARY_NAME if output == Path(".") else output, isdir=True)
        console.log(f"Writing to {fmt_filename(dirpath)}.")
        dirpath.mkdir()
        write_library(events, dataset, dirpath)

    console.print("\nDone.\n")
    # fmt: on
    return


DEFAULT_CONFIG_NAME = "mercury-config.yml"


@cli.command()
@click.argument(
    "output",
    type=click.Path(dir_okay=True, path_type=Path),
)
@click.pass_context
def drop(ctx: click.Context, output: Path):
    """Saves a yaml configuration default stub."""
    console = init_console(with_logo=False)
    config_text = (
        write_yaml(read_yaml(DEFAULT_CONFIG))  # removes comments
        if ctx.obj["quiet"]
        else DEFAULT_CONFIG
    )
    filepath = unused_path(
        output if not output.is_dir() else output / DEFAULT_CONFIG_NAME
    )
    with open(filepath, "w") as file:
        file.write(config_text)
    console.print(f"Created configuration file {fmt_filename(filepath)} :sparkles:.")


def sha1_hash(path: Path):
    with open(path, "rb") as f:
        # noinspection PyTypeChecker
        hsh = hashlib.file_digest(f, "sha1").hexdigest()
    return hsh


DEFAULT_EVENT_NAME = Path("event.fits")
DEFAULT_MERGE_NAME = ".mercury-merge.yaml"


@cli.command()
@click.argument(
    "input_directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.pass_context
def merge(ctx: click.Context, input_directory: Path):
    """Merge a library of results into its dataset."""
    from shutil import copy

    def fname(evtype: str, uuid_hex: str, uuid_substr: int = 4):
        index_uuid_substring = uuid_hex[:uuid_substr]
        return "".join(
            [
                DEFAULT_EVENT_NAME.stem,
                f"-{index_uuid_substring}",
                f"-{evtype}",
                DEFAULT_EVENT_NAME.suffix,
            ]
        )

    merge_path = input_directory / DEFAULT_MERGE_NAME
    # make sure the result directory was not already merged.
    if merge_path.is_file():
        raise click.UsageError("Input directory has already been merged.")

    # make sure an index file exists.
    index_path = input_directory / INDEX_FILENAME
    if not index_path.is_file():
        raise click.UsageError(
            f"Input directory does not contain a '{INDEX_FILENAME}' file."
        )
    with open(index_path, "r") as f:
        index = read_yaml(f)

    index_fmap = index["mappings"]
    files = [input_directory / file for file in index_fmap.keys()]
    roots = [Path(index_fmap[file]["root"]) for file in index_fmap.keys()]

    # checks that all result files and target directories exists
    if not all([f.exists() for f in [*files, *roots]]):
        raise click.FileError(
            f"Some of the indexed files or their root can not be found."
        )

    types = [index_fmap[file]["type"] for file in index_fmap.keys()]

    # move to destination folder and store info on moved files in a dict
    console = init_console(with_logo=False)
    merge_index = {}
    for file, root, event_type in zip(files, roots, types):
        # TODO: using `unused_path` we avoid potential clashes due to a repeated
        #  uuid substring but do not pad the filenames so the directory content
        #  can look messy if more than 10 files are generated. shall eventually
        #  find a better way to name these files.
        copy(file, dst := unused_path(root / fname(event_type, index["uuid"])))
        merge_index[str(file)] = {"dst": str(dst), "hash": sha1_hash(file)}

    # save infos to a yaml file that can be used to revert the merge
    with open(input_directory / DEFAULT_MERGE_NAME, "w") as f:
        write_yaml(merge_index, f)
    console.print(f"Merge complete :sparkles:.")


@cli.command()
@click.argument(
    "input_directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.pass_context
def clean(ctx: click.Context, input_directory: Path):
    """Undoes a merge, removing event data file from a dataset."""
    from os import remove

    # make sure the result directory was not already merged.
    merge_index_path = input_directory / DEFAULT_MERGE_NAME
    if not merge_index_path.is_file():
        raise click.UsageError(
            "The input directory is not a result folder or has not been merged yet."
        )

    with open(merge_index_path, "r") as f:
        merge_index = read_yaml(f)

    files = [Path(merge_index[file]["dst"]) for file in merge_index.keys()]
    hashes = [merge_index[file]["hash"] for file in merge_index.keys()]

    console = init_console(with_logo=False)
    for file, sha1hash in zip(files, hashes):
        if not file.is_file():
            warnings.warn(f"Skipping {file}. This file no longer exist.")
        elif sha1hash != sha1_hash(file):
            warnings.warn(f"Skipping {file}. Uncompatible hashes.")
        else:
            remove(file)

    if not any([file.is_file() for file in files]):
        console.print("All merged files removed. Cleaning complete :sparkles:!")
        remove(merge_index_path)
    else:
        console.print("Cleaning complete :sparkles:.")
