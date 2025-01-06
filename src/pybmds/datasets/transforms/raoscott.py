from enum import StrEnum
from typing import ClassVar, NamedTuple

import numpy as np
import pandas as pd

from ..dichotomous import DichotomousDataset


class Species(StrEnum):
    rat = "rat"
    mouse = "mouse"
    rabbit = "rabbit"


class Design(StrEnum):
    least_square = "Least Squares"
    orthogonal = "Orthogonal"


class P(NamedTuple):
    a: float
    b: float
    sigma: float


class RaoScott:
    adjustment_parameters: ClassVar = {
        (Species.mouse, Design.least_square): P(1.5938, 0.2866, 0.2078),
        (Species.mouse, Design.orthogonal): P(1.6943, 0.3132, 0.1863),
        (Species.rat, Design.least_square): P(1.6852, 0.331, 0.1248),
        (Species.rat, Design.orthogonal): P(1.8327, 0.369, 0.109),
        (Species.rabbit, Design.least_square): P(1.0582, 0.2397, 0.1452),
        (Species.rabbit, Design.orthogonal): P(1.1477, 0.2739, 0.1299),
    }

    def __init__(self, dataset: DichotomousDataset, species: Species):
        self.dataset = dataset
        self.species = species
        self.df = self.calculate()

    def calculate(self) -> pd.DataFrame:
        df = pd.DataFrame(
            {
                "dose": self.dataset.doses,
                "incidence": self.dataset.incidences,
                "n": self.dataset.ns,
            }
        )
        df["fraction_affected"] = df.incidence / df.n

        p = self.adjustment_parameters[(self.species, Design.least_square)]
        df["ls_design"] = np.exp(p.a + (p.b * np.log(df.fraction_affected)) + (0.5 * p.sigma))

        p = self.adjustment_parameters[(self.species, Design.orthogonal)]
        df["o_design"] = np.exp(p.a + (p.b * np.log(df.fraction_affected)) + (0.5 * p.sigma))

        df["avg_design"] = df[["ls_design", "o_design"]].mean(axis=1)
        df["scaled_incidence"] = df.incidence / df.avg_design
        df["scaled_n"] = df.n / df.avg_design
        return df

    def to_df(self) -> pd.DataFrame:
        return self.df
