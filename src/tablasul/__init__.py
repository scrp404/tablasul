"""Tablasul: First-principles PyWFS, Coronagraph, and Science PSF Data Simulator."""

from .core.field import ComplexElectricField
from .core.aperture import Aperture, CircularAperture, SubaruAperture
from .core.source import (
    Source,
    PointSource,
    CompanionSource,
    CalibrationSource,
    FilterBand,
)
from .core.atmosphere import KolmogorovScreen, FrozenFlowAtmosphere
from .optics.pywfs import PyramidWFS
from .optics.coronagraph import LyotCoronagraph
from .optics.ncpa import StaticPhaseNCPA, ZernikeNCPA
from .detector.camera import DetectorModel
from .engine.simulator import OpticalSystem, SimulationResult
from .engine.sync import SynchronizedCubeEngine, SynchronizedDataset
from .io.export import export_fits, export_npz

__version__ = "0.1.0"
__all__ = [
    "ComplexElectricField",
    "Aperture",
    "CircularAperture",
    "SubaruAperture",
    "Source",
    "PointSource",
    "CompanionSource",
    "CalibrationSource",
    "FilterBand",
    "KolmogorovScreen",
    "FrozenFlowAtmosphere",
    "PyramidWFS",
    "LyotCoronagraph",
    "StaticPhaseNCPA",
    "ZernikeNCPA",
    "DetectorModel",
    "OpticalSystem",
    "SimulationResult",
    "SynchronizedCubeEngine",
    "SynchronizedDataset",
    "export_fits",
    "export_npz",
]
