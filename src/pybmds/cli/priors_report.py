import argparse
import os
import sys
from datetime import datetime
from io import StringIO
from pathlib import Path

from ..constants import DistType, PriorClass
from ..datasets import ContinuousDataset, DichotomousDataset
from ..models import continuous, dichotomous


def write_model(f: StringIO, ModelClass: type[dichotomous.BmdModel]):
    f.write(f"### {ModelClass.__name__}\n\n")


def write_break(f: StringIO):
    f.write(f'{"-"*80}\n\n')


def write_settings(f: StringIO, model: dichotomous.BmdModel, settings: dict):
    f.write("\n".join(f"* {k}: {v!r}" for k, v in settings.items()) + "\n\n")
    f.write(model.settings.priors.report_tbl() + "\n\n")


def dichotomous_priors(f: StringIO):
    dichotomous_dataset = DichotomousDataset(
        doses=[0, 1.96, 5.69, 29.75], ns=[75, 49, 50, 49], incidences=[5, 1, 3, 14]
    )

    def _print_d_model(ModelClass: type[dichotomous.BmdModelDichotomous], restricted: bool):
        write_model(f, ModelClass)

        # print unrestricted
        settings = {"priors": PriorClass.frequentist_unrestricted}
        model = ModelClass(dataset=dichotomous_dataset, settings=settings)
        write_settings(f, model, settings)

        # print restricted
        if restricted:
            settings = {"priors": PriorClass.frequentist_restricted}
            model = ModelClass(dataset=dichotomous_dataset, settings=settings)
            write_settings(f, model, settings)

        # print bayesian
        settings = {"priors": PriorClass.bayesian}
        model = ModelClass(dataset=dichotomous_dataset, settings=settings)
        write_settings(f, model, settings)

        write_break(f)

    f.write("## Dichotomous\n\n")
    _print_d_model(dichotomous.LogLogistic, True)
    _print_d_model(dichotomous.Gamma, True)
    _print_d_model(dichotomous.Logistic, False)
    _print_d_model(dichotomous.Probit, False)
    _print_d_model(dichotomous.QuantalLinear, False)
    _print_d_model(dichotomous.LogProbit, True)
    _print_d_model(dichotomous.Weibull, True)
    _print_d_model(dichotomous.Multistage, True)
    _print_d_model(dichotomous.DichotomousHill, True)


def continuous_priors(f: StringIO):
    continuous_dataset = ContinuousDataset(
        doses=[0, 10, 50, 150, 400],
        ns=[10, 10, 10, 10, 10],
        means=[10, 20, 30, 40, 50],
        stdevs=[1, 2, 3, 4, 5],
    )

    def print_c_model(ModelClass: type[continuous.BmdModelContinuous], settings: dict):
        model = ModelClass(dataset=continuous_dataset, settings=settings)
        write_settings(f, model, settings)

    f.write("## Continuous\n\n")

    write_model(f, continuous.ExponentialM3)
    for settings in [
        dict(priors=PriorClass.frequentist_restricted, is_increasing=True),
        dict(priors=PriorClass.frequentist_restricted, is_increasing=False),
        dict(priors=PriorClass.bayesian),
    ]:
        print_c_model(continuous.ExponentialM3, settings)
    write_break(f)

    write_model(f, continuous.ExponentialM5)
    for settings in [
        dict(priors=PriorClass.frequentist_restricted, is_increasing=True),
        dict(priors=PriorClass.frequentist_restricted, is_increasing=False),
        dict(priors=PriorClass.bayesian),
    ]:
        print_c_model(continuous.ExponentialM5, settings)
    write_break(f)

    write_model(f, continuous.Hill)
    print_c_model(continuous.Hill, dict(priors=PriorClass.frequentist_restricted))
    print_c_model(continuous.Hill, dict(priors=PriorClass.frequentist_unrestricted))
    print_c_model(continuous.Hill, dict(priors=PriorClass.bayesian))
    write_break(f)

    write_model(f, continuous.Linear)
    # fmt: off
    for settings in [
        dict(priors=PriorClass.frequentist_unrestricted, disttype=DistType.normal),
        dict(priors=PriorClass.frequentist_unrestricted, disttype=DistType.normal_ncv),

        dict(priors=PriorClass.bayesian),
    ]:  # fmt: on
        print_c_model(continuous.Linear, settings)
    write_break(f)

    write_model(f, continuous.Polynomial)
    # fmt: off
    for settings in [
        dict(priors=PriorClass.frequentist_unrestricted, disttype=DistType.normal),
        dict(priors=PriorClass.frequentist_unrestricted, disttype=DistType.normal_ncv),

        dict(priors=PriorClass.frequentist_restricted, disttype=DistType.normal, is_increasing=True),
        dict(priors=PriorClass.frequentist_restricted, disttype=DistType.normal, is_increasing=False),
        dict(priors=PriorClass.frequentist_restricted, disttype=DistType.normal_ncv, is_increasing=True),
        dict(priors=PriorClass.frequentist_restricted, disttype=DistType.normal_ncv, is_increasing=False),

        dict(priors=PriorClass.bayesian),
    ]:  # fmt: on
        print_c_model(continuous.Polynomial, settings)
    write_break(f)

    write_model(f, continuous.Power)
    # fmt: off
    for settings in [
        dict(priors=PriorClass.frequentist_unrestricted, disttype=DistType.normal),
        dict(priors=PriorClass.frequentist_unrestricted, disttype=DistType.normal_ncv),

        dict(priors=PriorClass.frequentist_restricted, disttype=DistType.normal),
        dict(priors=PriorClass.frequentist_restricted, disttype=DistType.normal_ncv),

        dict(priors=PriorClass.bayesian),
    ]:  # fmt: on
        print_c_model(continuous.Power, settings)
    write_break(f)


def create_report() -> StringIO:
    f = StringIO()
    f.write(f"# BMDS priors report:\n\nGenerated on: {datetime.now()}\n\n")
    dichotomous_priors(f)
    continuous_priors(f)
    return f


def main():
    parser = argparse.ArgumentParser(description="Generate a report on BMDS priors settings")
    parser.add_argument(
        "filename", nargs="?", help="Optional; output file. If empty, writes to stdout."
    )
    args = parser.parse_args()
    report = create_report()
    if args.filename:
        path = (Path(os.curdir) / sys.argv[1]).expanduser().resolve().absolute()
        sys.stdout.write(f"Writing output to: {path}")
        path.write_text(report.getvalue())
    else:
        sys.stdout.write(report.getvalue())
