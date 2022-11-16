import typer

from . import work

app = typer.Typer()


@app.command()
def hello(name: str = "world"):
    typer.secho(f"Hello {name}!", fg="red")


@app.command()
def bottles(num: int = 10, beverage: str = "coke"):
    work.bottles(num, beverage)
