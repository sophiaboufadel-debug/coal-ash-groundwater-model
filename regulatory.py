"""
regulatory.py
=============
Checklist section 38: kept deliberately separate from the transport physics
so changing a regulatory threshold never changes a simulated concentration
-- only the exceedance ratio computed from it.
"""
from dataclasses import dataclass

from reference_data import REGULATORY_STANDARDS, RegulatoryStandard


@dataclass
class ExceedanceResult:
    contaminant: str
    predicted_mg_l: float
    standard_mg_l: float
    standard_kind: str
    ratio: float
    exceeds: bool
    note: str


def exceedance_ratio(contaminant: str, predicted_mg_l: float) -> ExceedanceResult:
    """ER = C_predicted / C_limit   (checklist #38)."""
    if contaminant not in REGULATORY_STANDARDS:
        raise KeyError(f"No regulatory standard on file for '{contaminant}'. "
                        f"Add it to reference_data.REGULATORY_STANDARDS.")
    std: RegulatoryStandard = REGULATORY_STANDARDS[contaminant]
    ratio = predicted_mg_l / std.value_mg_l if std.value_mg_l > 0 else float("inf")
    return ExceedanceResult(
        contaminant=contaminant,
        predicted_mg_l=predicted_mg_l,
        standard_mg_l=std.value_mg_l,
        standard_kind=std.kind,
        ratio=ratio,
        exceeds=ratio > 1.0,
        note=std.note,
    )


def exceedance_table(predicted_by_contaminant: dict) -> list:
    """Batch version: {contaminant: predicted_mg_l} -> list[ExceedanceResult]."""
    return [exceedance_ratio(c, v) for c, v in predicted_by_contaminant.items() if c in REGULATORY_STANDARDS]
