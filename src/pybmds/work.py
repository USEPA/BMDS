import time
from random import random

import typer

from .bleep import add


def bottles(num: int, beverage: str):
    for i in range(num, 0, -1):
        typer.secho(
            f"{i} bottles of {beverage} on the wall, {i} bottles of {beverage}, take one down...",
            fg="green",
        )
        time.sleep(random() / 3)
    typer.secho(f"No more bottles of {beverage} on the wall!", fg="green")


def super_add(a: int, b: int) -> int:
    # a silly function for an example unit-test
    return a + b


def double_add(a: int, b: int) -> int:
    return add(a, b) * 2
