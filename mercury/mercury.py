from pathlib import Path
from typing import Callable, Sequence

import click
from rich.console import Console
from schema import And  # type: ignore[import-untyped]
from schema import Optional  # type: ignore[import-untyped]
from schema import Schema  # type: ignore[import-untyped]
from schema import SchemaError  # type: ignore[import-untyped]
from schema import Use  # type: ignore[import-untyped]
from yaml import YAMLError
from yaml import dump as write_yaml
from yaml import safe_load as read_yaml

import hbstools as hbs

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
 # |               | Python PF+DES | Python BFT | C PF+SES | C BFT |
 # |---------------|---------------|------------|----------|-------|
 # | threshold_std |       ✓       |      ✓     |     ✓    |   ✓   |
 # | mu_min        |       ✓       |      ✓     |     ✓    |   ✓   |
 # | alpha         |       ✓       |      ✓     |     ✓    |   ✓   |
 # | beta          |       ✓       |      ✓     |          |       |
 # | m             |       ✓       |      ✓     |     ✓    |   ✓   |
 # | t_max         |       ✓       |      ✓     |          |       |
 # | sleep         |       ✓       |      ✓     |     ✓    |   ✓   |
 # | majority      |               |      ✓     |          |   ✓   |
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


def _default_config() -> tuple[Callable, Callable]:
    def show() -> str:
        """Returns default config"""
        if called_yet[0]:
            return DEFAULT_CONFIG
        else:
            called_yet[0] = True
            return DEFAULT_CONFIG

    def tell() -> bool:
        """A sentinel. It tells if default config was ever called"""
        return called_yet[0]

    called_yet = [False]
    return show, tell


default_config, used_default_config = _default_config()


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


def parse_user_else_default_config(config_path: Path | None) -> dict:
    """Opens the configuration file, if provided, or get values from default."""
    if config_path is None:
        configuration = read_yaml(default_config())
    else:
        with open(config_path, "r") as stream:
            configuration = read_yaml(stream)
    return configuration


def crawler(
    directory: Path | str,
    targets: Sequence[str],
    recursion_limit: int = 1,
) -> set[Path]:
    """Crawls through folder looking for subfolders containing a specific set of files.
    if `recursion_limit=0`  we only check the present folder, ignoring its subdirectories.
    We return a set, which has no meaningful order."""

    def get_subdirs(d: Path) -> list[Path]:
        return [f for f in d.iterdir() if f.is_dir()]

    def directory_contains(d: Path, ts) -> bool:
        for file in ts:
            if not d.joinpath(file).is_file():
                return False
        return True

    def check_subdirs(
        subdirs: list[Path], ts, rlim: int, recursion_acc: int, acc: set[Path]
    ):
        if not subdirs:
            return []
        car, *cdr = subdirs
        check_dir(car, ts, rlim, recursion_acc + 1, acc)
        check_subdirs(cdr, ts, rlim, recursion_acc, acc)

    def check_dir(d: Path, ts, rlim, recursion_acc: int, acc: set[Path]) -> set[Path]:
        if recursion_acc == rlim + 1:
            return acc
        if directory_contains(d, ts):
            acc.add(d)
        check_subdirs(get_subdirs(d), ts, rlim, recursion_acc, acc)
        return acc

    return check_dir(Path(directory), targets, recursion_limit, 0, set())


def fmt_filename(filename: str | Path) -> str:
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
    """Validates user configuration option and records it."""
    # checks that the YAML config is well written, with well-defined values.
    try:
        config = config_schema.validate(parse_user_else_default_config(config_path))
    except YAMLError:
        raise click.BadParameter("Cannot parse YAML configuration file.")
    except SchemaError as error:
        raise click.BadParameter(f"Wrong inputs in config file: \n - {error}")

    # check if we have an algorithm matching the user configuration.
    algorithm_params = config["algorithm_params"]
    try:
        algorithm_class = hbs.trigger.trigger_match(algorithm_params)
    except ValueError:
        raise click.BadParameter("Cannot find an algorithm matching the configuration.")

    # store the algorithm name so that we can show it during execution.
    ctx.obj["search_algoname"] = str(algorithm_class(**algorithm_params))
    ctx.obj["search_config"] = None if used_default_config() else config_path
    return config


def unused_path(file: Path, num: int = 1, isdir: bool = False) -> Path:
    """Return an unused path with same stem prefix as `file` and incremental suffix.
    if `num` is set to 1, tested alternatives are: `path.fits`, `path-1.fits`,
    `path-2.fits` and so on.
    if `num` is set to 0, tested alternatives are `path.fits`, `path-0.fits`,
    `path-1.fits`..`"""
    if (not file.is_file() and not isdir) or (not file.is_dir() and isdir):
        return file
    parts = file.stem.split("-")
    *tail, head = parts
    stem = "-".join(tail) if head.isdigit() else "-".join(parts)
    return unused_path(Path(file).parent.joinpath(f"{stem}-{num}{file.suffix}"), num + 1, isdir)


DEFAULT_CATALOG_NAME = "mercury-results.fits"
DEFAULT_LIBRARY_NAME = "mercury-results/"


@cli.command()
@click.argument(
    "input_dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
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
    "In `library` mode, a folder or fits file is generated"
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
    input_dir: Path,
    output: Path,
    configuration: dict,
    mode: str,
    recursion_limit: int,
):
    """Searches transients from data in a folder.
    The algorithm used for this search is called Poisson-FOCuS (Ward, 2022; Dilillo, 2024).
    The search is configurable using a YAML configuration, see mercury's 'drop' command.
    """
    # fmt: off
    console = init_console(with_logo=True)
    console.print(f"Welcome, this is [bold]mercury.search[/].\n")
    config_file = "default" if (f := ctx.obj["search_config"]) is None else f
    console.log(f"Loaded {fmt_filename(config_file)} configuration.")
    console.log(f"Algorithm {ctx.obj['search_algoname']} matches the configuration.")

    search_targets = ["gti.fits", "out_s_cl.evt", "out_x_cl.evt"]
    data_paths = crawler(input_dir, search_targets, recursion_limit)
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
        hbs.io.write_catalog(events, filepath)
    elif mode == "library":
        dirpath = unused_path(output / DEFAULT_LIBRARY_NAME if output == Path(".") else output, isdir=True)
        dirpath.mkdir()
        console.log(f"Writing to {fmt_filename(dirpath)}.")
        hbs.io.write_library(events, dataset, dirpath)

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
        write_yaml(read_yaml(default_config()))
        if ctx.obj["quiet"]
        else default_config()
    )
    filepath = unused_path(output if not output.is_dir() else output / DEFAULT_CONFIG_NAME)
    with open(filepath, "w") as file:
        file.write(config_text)
    console.print(f"Created configuration file '{fmt_filename(filepath)}' :sparkles:.")


@cli.command()
@click.argument(
    "input_directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.pass_context
def merge(ctx: click.Context, input_directory: Path):
    """Merge a library of results into a dataset."""
    from shutil import copy

    console = init_console(with_logo=False)
    index_path = input_directory / hbs.io.INDEX_FILENAME
    if not index_path.is_file():
        raise click.UsageError(
            f"The input folders does not contain a '{hbs.io.INDEX_FILENAME}' file."
        )
    with open(index_path, "r") as f:
        index = read_yaml(f)

    files, roots = zip(*index.items())

    assert len(files) == len(roots)
    if not all([f.exists() for f in map(Path, [*files, *roots])]):
        raise click.FileError(
            f"Some of the indexed files or their root can not be found."
        )
    if not all([not (Path(root) / f).exists() for f, root in zip(files, roots)]):
        print([(Path(root) / f).exists() for f, root in zip(files, roots)])
        print([Path(root) / f for f, root in zip(files, roots)])
        raise click.FileError(
            f"The files were already merged, or the roots folder contain event files."
        )

    for file, root in zip(files, roots):
        copy(file, root)
        console.print(f"Copied {fmt_filename(Path(file).stem)} to {fmt_filename(Path(root))}.")

    console.print(f"\nMerge complete :sparkles:.")
