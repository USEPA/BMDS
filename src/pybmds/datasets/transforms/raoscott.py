from enum import StrEnum

import numpy as np
import pandas as pd

from ..dichotomous import DichotomousDataset


class Species(StrEnum):
    rat = "Rat"
    mouse = "Mouse"
    rabbit = "Rabbit"


class Design(StrEnum):
    least_square = "Least Squares"
    orthogonal = "Orthogonal"


params = {
    (Species.mouse, Design.least_square): (1.5938, 0.2866, 0.2078),
    (Species.mouse, Design.orthogonal): (1.6943, 0.3132, 0.1863),
    (Species.rat, Design.least_square): (1.6852, 0.331, 0.1248),
    (Species.rat, Design.orthogonal): (1.8327, 0.369, 0.109),
    (Species.rabbit, Design.least_square): (1.0582, 0.2397, 0.1452),
    (Species.rabbit, Design.orthogonal): (1.1477, 0.2739, 0.1299),
}


def rao_scott(dataset: DichotomousDataset, species: Species) -> pd.DataFrame:
    df = pd.DataFrame(
        {
            "dose": dataset.doses,
            "incidence": dataset.incidences,
            "n": dataset.ns,
        }
    )
    df["fraction_affected"] = df.incidence / df.n
    a, b, sigma = params[(species, Design.least_square)]
    df["ls_design"] = np.exp(a + (b * np.log(df.fraction_affected)) + (0.5 * sigma))
    a, b, sigma = params[(species, Design.orthogonal)]
    df["o_design"] = np.exp(a + (b * np.log(df.fraction_affected)) + (0.5 * sigma))
    df["avg_design"] = df[["ls_design", "o_design"]].mean(axis=1)
    df["scaled_incidence"] = df.incidence / df.avg_design
    df["scaled_n"] = df.n / df.avg_design
    return df
