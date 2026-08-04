"""
flow.py
=======
Checklist sections 2-9, 22-25: water movement. Darcy's law, pore velocity,
travel time, and the vertical-infiltration / horizontal-transport split
(contaminant first moves down through the ash and underlying soil beneath
the pond, then laterally once it reaches the water table).
"""
from dataclasses import dataclass
import numpy as np

from units import SECONDS_PER_YEAR


def darcy_flux(k_m_s: float, gradient: float) -> float:
    """q = K * i   (checklist #5). Returns flux in m/s."""
    return k_m_s * gradient


def darcy_flow_rate(k_m_s: float, area_m2: float, gradient: float) -> float:
    """Q = K * A * i   (checklist #5). Returns m^3/s."""
    return k_m_s * area_m2 * gradient


def pore_velocity(darcy_flux_m_s: float, effective_porosity: float) -> float:
    """v = q / n_e   (checklist #6). Returns m/s."""
    if effective_porosity <= 0:
        raise ValueError("effective_porosity must be positive")
    return darcy_flux_m_s / effective_porosity


def hydraulic_gradient(h1_m: float, h2_m: float, length_m: float) -> float:
    """i = (h1 - h2) / L   (checklist #4)."""
    if length_m <= 0:
        raise ValueError("length_m must be positive")
    return (h1_m - h2_m) / length_m


def unretarded_travel_time_years(distance_m: float, velocity_m_s: float) -> float:
    """t = x / v   (checklist #7). velocity_m_s must be > 0."""
    if velocity_m_s <= 0:
        raise ValueError("velocity_m_s must be positive")
    return (distance_m / velocity_m_s) / SECONDS_PER_YEAR


@dataclass
class VerticalInfiltration:
    """Checklist #22-23: water and contaminant moving DOWN beneath the pond,
    through the ash and then the underlying soil, before reaching the water
    table. Modeled as two conductivities in series (ash layer, then soil
    layer) using a harmonic-mean effective vertical conductivity -- the
    standard treatment for vertical flow across layered media
    (Freeze & Cherry, Groundwater, 1979, Ch. 2)."""
    k_ash_m_s: float
    thickness_ash_m: float
    k_soil_m_s: float
    thickness_soil_m: float

    @property
    def effective_k_v_m_s(self) -> float:
        total_thickness = self.thickness_ash_m + self.thickness_soil_m
        return total_thickness / (
            self.thickness_ash_m / self.k_ash_m_s + self.thickness_soil_m / self.k_soil_m_s
        )

    def flux_m_s(self, head_difference_m: float) -> float:
        total_thickness = self.thickness_ash_m + self.thickness_soil_m
        gradient = head_difference_m / total_thickness
        return darcy_flux(self.effective_k_v_m_s, gradient)

    def infiltration_volume_m3_per_year(self, head_difference_m: float, pond_area_m2: float) -> float:
        flux = self.flux_m_s(head_difference_m)
        return max(flux, 0.0) * pond_area_m2 * SECONDS_PER_YEAR


@dataclass
class HorizontalTransport:
    """Checklist #24-25: once in the aquifer, contaminant moves laterally."""
    k_h_m_s: float
    gradient: float
    effective_porosity: float

    @property
    def darcy_flux_m_s(self) -> float:
        return darcy_flux(self.k_h_m_s, self.gradient)

    @property
    def pore_velocity_m_s(self) -> float:
        return pore_velocity(self.darcy_flux_m_s, self.effective_porosity)

    @property
    def pore_velocity_m_year(self) -> float:
        return self.pore_velocity_m_s * SECONDS_PER_YEAR

    def flow_rate_m3_s(self, cross_sectional_area_m2: float) -> float:
        return darcy_flow_rate(self.k_h_m_s, cross_sectional_area_m2, self.gradient)
