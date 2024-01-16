from math import log10
from pathlib import Path
from typing import Callable

import click
from rich.console import Console
from schema import And  # type: ignore[import-untyped]
from schema import Schema  # type: ignore[import-untyped]
from schema import SchemaError  # type: ignore[import-untyped]
from schema import Use  # type: ignore[import-untyped]
from yaml import dump as write_yaml
from yaml import safe_load as read_yaml
from yaml import YAMLError

from hbstools.io import write_ttis_to_fits
from hbstools.search import Search

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
 # The `threshold` keyword sets the confidence level for detection.
 # The higher the threshold, the less false positive but the more the false negatives.
 # The threshold is expressed in units of standard deviations.
 # Must be greater than 0.
 threshold: 5.5

 # The `mu_min` key set the focus parameter for killing old changepoints 
 # which most likely will never result in a trigger. Keep it below 1.5.
 # Must be greater or equal than 1.0. Disabled if it equals 1.0.
 mu_min: 1.1
 
 # The `alpha` keyword sets the background smoothing characteristic time \tau, 
 # where \tau = (binning/0.005).
 # Must be greater than 0.
 alpha: 0.005

 # The `beta` keyword sets the trend component for background estimate.
 # The algorithm may become unstable when `beta` is set: leave it `0.0` unless
 # you have good reasons to use it.
 # Must be non-negative.
 beta: 0.0

 # The `m` keyword will prevent the most recent observed counts to be used for 
 # background estimation. This prevents background estimate to be "polluted"
 # by real transients. 
 # It is expressed in units of bin-steps. This means that if `binning` is set
 # to 0.1 and `m` is set to 40, the algorithm won't use the latest 4.0 s
 # of data for background estimate. 
 # Must be a non-negative integer.
 m: 40

 # The `t_max` parameter tells the algorithm to kill old changepoint.
 # It is a good idea to keep it equal to `m`. 
 # It is expressed as a bin-step, see `m` or `skip`.
 # Must be an integer greater than 0.
 t_max: 40

 # The algorithms stays idle for a while before starting its operations.
 # This help forming a good estimate of the background. 
 # The `sleep` parameters sets for how long this idle period lasts.
 # The `sleep` parameter is expressed as a bin-step, see `m` or `skip`.
 # Must be an integer greater than `m`.
 sleep: 1600

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
        """Tells if default config was ever called"""
        return called_yet[0]

    called_yet = [False]
    return show, tell


default_config, used_default_config = _default_config()


class ConfigSchema(Schema):
    """A schema used in validation of mercury's configuration files."""

    def validate(self, data, _is_event_schema=True):
        data = super(ConfigSchema, self).validate(data, _is_event_schema=False)
        if (
            _is_event_schema
            and data["algorithm_params"]["sleep"] <= data["algorithm_params"]["m"]
        ):
            raise SchemaError("`sleep` must be strictly greater than `m`.")
        return data


config_schema = ConfigSchema(
    {
        "binning": And(
            Use(float),
            lambda b: b > 0,
            error="`binning` must be greater than zero.",
        ),
        "energy_lims": And(
            Use(tuple[float, float]),
            lambda l: 0 <= l[0] < l[1],
        ),
        "skip": And(
            lambda s: isinstance(s, int) & (s >= 0),
            error="`skip` must be a non-negative integer",
        ),
        "algorithm_params": {
            "threshold": And(
                Use(float),
                lambda t: t > 0,
                error="`threshold` must be positive.",
            ),
            "alpha": And(
                Use(float),
                lambda a: a > 0,
                error="`alpha` must be positive",
            ),
            "beta": And(
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
                lambda m: isinstance(m, int) & (m >= 0),
                error="`m` must be a non negative integer.",
            ),
            "t_max": And(
                lambda t: isinstance(t, int) & (t > 0),
                error="`t_max must be an integer greater than zero",
            ),
            "sleep": And(
                lambda t: isinstance(t, int) & (t >= 0),
                error="`sleep` must be non negative",
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
    target: list[str],
    recursion_limit: int = 2,
) -> list[Path]:
    """Crawls through folder looking for subfolders containing a specific set of files."""

    def get_subdirs(dir_: Path) -> list[Path]:
        return [d for d in dir_.iterdir() if d.is_dir()]

    def directory_contains(dir_: Path, target: list[Path]) -> bool:
        for file in target:
            if not dir_.joinpath(file).is_file():
                return False
        return True

    def check_subdirs(subdirs, lf, reclim, recnum, acc):
        if not subdirs:
            return []
        car, *cdr = subdirs
        check_dir(car, lf, reclim, recnum + 1, acc)
        check_subdirs(cdr, lf, reclim, recnum, acc)

    def check_dir(dir_, lf, reclim, recnum: int, acc: list[Path]) -> list[Path]:
        if directory_contains(dir_, lf):
            acc.append(dir_)
        if recnum == reclim:
            return acc
        check_subdirs(get_subdirs(dir_), lf, reclim, recnum, acc)
        return acc

    return check_dir(Path(directory), target, recursion_limit, 0, [])


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


def search_record_input(ctx: click.Context, param: click.Option, value: Path) -> Path:
    """Records user input argument"""
    ctx.obj["search_input"] = value
    return value


def search_validate_config(
    ctx: click.Context, param: click.Option, value: Path | None
) -> dict:
    """Validates user configuration option and records it."""

    def validate_config(config_path: Path | None) -> dict:
        config = parse_user_else_default_config(config_path)
        config_schema.validate(config)
        return config

    try:
        configuration = validate_config(value)
    except YAMLError:
        raise click.BadParameter("Cannot parse YAML configuration file.")
    except SchemaError as error:
        raise click.BadParameter(f"Wrong inputs in config file:\n - {error}")

    ctx.obj["search_config"] = None if used_default_config() else value
    return configuration


def validate_output(output: Path, filename: str) -> Path:
    """Checks output choice and finds alternatives when appropriate."""

    def find_valid_name(file, num):
        if not file.exists():
            return file
        else:
            next_name = f"{file.stem[:-(int(log10(num)) + 2)] if num > 1 else file.stem}-{num}{file.suffix}"
            return find_valid_name(Path(file).parent.joinpath(next_name), num + 1)

    if output.is_dir():
        default_file = Path(output).joinpath(filename)
        if default_file.is_file():
            return find_valid_name(default_file, 1)
        else:
            return default_file
    elif output.is_file():
        raise FileExistsError()
    else:
        return output


def search_validate_output(
    ctx: click.Context, param: click.Option, value: Path
) -> Path:
    """Validate user output choice and records it"""
    try:
        output = validate_output(
            value if value is not None else ctx.obj["search_input"],
            "mercury-results.fits",
        )
    except FileExistsError:
        raise click.BadParameter(
            f"Can not create file, a file '{value}' already exists."
        )
    ctx.obj["seach_output"] = value
    return output


@cli.command()
@click.argument(
    "input_directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    callback=search_record_input,
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
    "-o",
    "--output",
    "output_path",
    type=click.Path(dir_okay=True, path_type=Path),
    callback=search_validate_output,
    help="Path to a file or a directory where results are saved. "
    "If not provided, the output are saved in the input directory. ",
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
    input_directory: Path,
    configuration: dict,
    output_path: Path,
    recursion_limit: int,
):
    """Searches transients from data in a folder.
    The algorithm used for this search is called Poisson-FOCuS (Ward, 2022; Dilillo, 2024).
    The search is configurable using a YAML configuration, see mercury's 'drop' command.
    """

    def fmt_filename(filename: str | Path) -> str:
        return f"'[b]{filename}[/]'"

    console = init_console(with_logo=True)
    console.print(f"Welcome, this is [bold]mercury.search[/].\n")
    config_file = "default" if (f := ctx.obj["search_config"]) is None else f
    console.log(f"Loaded {fmt_filename(config_file)} configuration.")

    search_targets = ["gti.fits", "out_s_cl.evt", "out_x_cl.evt"]
    dataset = crawler(input_directory, search_targets, recursion_limit)
    if not dataset:
        console.print("\nFound no data. Exiting.\n")
        return
    else:
        console.log(f"Found {len(dataset)} dataset.")

    _search = Search(**configuration, console=console)
    ttis = _search(dataset)

    write_ttis_to_fits(ttis, output_path)
    console.log(f"Writing to {fmt_filename(output_path)}.")

    console.print("\nDone.\n")
    return


def drop_validate_output(ctx: click.Context, param: click.Option, value: Path) -> Path:
    """Validate user output choice and records it."""
    try:
        output = validate_output(value, "mercury-config.yml")
    except FileExistsError:
        raise click.BadParameter(
            f"Can not create file, a file '{value}' already exists."
        )
    ctx.obj["drop_output"] = value
    return output


@cli.command()
@click.argument(
    "output",
    type=click.Path(dir_okay=True, path_type=Path),
    callback=drop_validate_output,
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
    with open(output, "w") as file:
        file.write(config_text)
    console.print(f"[bold]Created configuration file '{output}' :sparkles:.")
