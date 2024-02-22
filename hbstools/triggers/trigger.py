from hbstools.triggers.poissonfocusdes import PoissonFocusDES
from hbstools.triggers.bft import Bft
from hbstools.triggers.poissonfocusses_cwrap import PoissonFocusSES_C
from hbstools.triggers.bft_cwrap import Bft_C


class TriggerAlgorithm:
    def __init__(
            self,
            algorithm: PoissonFocusDES | Bft | PoissonFocusSES_C | Bft_C,
            energy_lims: tuple,
            binning: float,
            algorithm_params: dict,
    ):
        self.algorithm = algorithm
        self.energy_lims = energy_lims
        self.binning = binning
        self.algorithm_params = algorithm_params

