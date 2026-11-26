"""Separate FITC-folate equilibrium component, not a CD19/CD22 construct model."""

from dataclasses import dataclass
from importlib.resources import files
import json
import math

import numpy as np

from .validation import number


@dataclass(frozen=True)
class AdapterBinding:
    """All concentrations nM; free adapter is an externally maintained reservoir."""

    car_kd: float = 0.03
    target_kd: float = 1.0
    car_total: float = 0.01
    target_total: float = 0.01

    def __post_init__(self):
        for name, value in vars(self).items():
            number(value, name, positive=True)

    @classmethod
    def published_affinities(cls):
        data = json.loads(
            files("cart_dynamics").joinpath("data/sources.json").read_text()
        )["lee2019"]
        return cls(
            car_kd=data["car_cell_binding_kd_nM"],
            target_kd=data["folate_receptor_binding_kd_nM"],
        )

    @property
    def peak_free_adapter(self):
        return math.sqrt(self.car_kd) * math.sqrt(self.target_kd)

    def species(self, adapter):
        number(adapter, "free adapter")
        # c = A / ((Kd_R+A)(Kd_F+A)); stable form avoids A^2 overflow.
        c = (adapter / (self.car_kd + adapter)) / (self.target_kd + adapter)
        x = c * (self.car_total + self.target_total)
        discriminant = 1 + 2 * x + (c * (self.car_total - self.target_total)) ** 2
        bridge = (2 * c * self.car_total * self.target_total) / (
            1 + x + math.sqrt(discriminant)
        )
        free_car = (self.car_total - bridge) / (1 + adapter / self.car_kd)
        free_target = (self.target_total - bridge) / (1 + adapter / self.target_kd)
        return {
            "bridge": bridge,
            "free_car": free_car,
            "free_target": free_target,
            "car_adapter": free_car * adapter / self.car_kd,
            "target_adapter": free_target * adapter / self.target_kd,
        }

    def bridge_fraction(self, adapter):
        return self.species(adapter)["bridge"] / min(self.car_total, self.target_total)


def free_adapter_profile(minutes, pulses, half_life):
    """Specified free-concentration jumps, not administered doses or a full PK model."""
    number(half_life, "half_life", positive=True)
    minutes = np.asarray(minutes, dtype=float)
    if not np.isfinite(minutes).all() or np.any(minutes < 0):
        raise ValueError("minutes must be finite and nonnegative")
    concentration = np.zeros_like(minutes)
    for time, amount in pulses:
        number(time, "pulse time")
        number(amount, "free-concentration jump")
        elapsed = np.maximum(minutes - time, 0)
        concentration += (
            amount * np.exp(-math.log(2) * elapsed / half_life) * (minutes >= time)
        )
    return concentration
