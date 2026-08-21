from types import SimpleNamespace

import pandas as pd

from pybmds.constants import BMDS_BLANK_VALUE
from pybmds.reporting.styling import (
    Report,
    _ma_model_bmd_triplet,
    df_to_table,
    parameter_summary_formatter,
    write_bayesian_table,
    write_cell,
    write_dataset_table,
)


def test_write_dataset_table(cdataset, cidataset, ddataset, nd_dataset):
    report = Report.build_default()
    write_dataset_table(report, cdataset, True)
    write_dataset_table(report, cdataset, False)
    write_dataset_table(report, cidataset, True)
    write_dataset_table(report, cidataset, False)
    write_dataset_table(report, ddataset, True)
    write_dataset_table(report, ddataset, False)
    write_dataset_table(report, nd_dataset, False)


def test_df_to_table_writes_dataframe_footnotes():
    report = Report.build_default()
    df = pd.DataFrame({"Model": ["A"], "BMD": [1.0]})
    df.attrs["footnotes"] = [
        "R-hat statistic is calculated only when more than 1 Markov chain is used."
    ]

    df_to_table(report, df)

    assert report.document.paragraphs[-1].text == (
        "R-hat statistic is calculated only when more than 1 Markov chain is used."
    )


def test_df_to_table_accepts_custom_formatter():
    report = Report.build_default()
    df = pd.DataFrame({"Parameter": ["alpha"], "Median": [0.000584]})

    df_to_table(report, df, formatter=parameter_summary_formatter)

    assert report.document.tables[0].cell(1, 1).text == "0.000584"


def test_write_cell_formats_blank_value():
    report = Report.build_default()
    table = report.document.add_table(1, 1)

    write_cell(table.cell(0, 0), BMDS_BLANK_VALUE, report.styles.tbl_body)

    assert table.cell(0, 0).text == "-"


def test_parameter_summary_formatter_keeps_small_values():
    assert parameter_summary_formatter("already formatted") == "already formatted"
    assert parameter_summary_formatter(float("nan")) == "-"
    assert parameter_summary_formatter(0.000584) == "0.000584"
    assert parameter_summary_formatter(-0.000584) == "-0.000584"
    assert parameter_summary_formatter(1.25) == "1.25"


def _bayesian_summary_headers(is_loud: bool) -> list[str]:
    report = Report.build_default()
    model = SimpleNamespace(
        results=SimpleNamespace(
            bmdl=1.0,
            bmd=2.0,
            bmdu=3.0,
            fit=SimpleNamespace(bic_equiv=4.0),
            gof=SimpleNamespace(residual=[0.1], roi=0.2),
        ),
        name=lambda: "Model 1",
    )
    ma = SimpleNamespace(
        has_results=True,
        models=[model],
        results=SimpleNamespace(
            priors=[0.5],
            posteriors=[0.5],
            bmdl=1.0,
            bmd=2.0,
            bmdu=3.0,
        ),
    )
    session = SimpleNamespace(
        models=[model],
        model_average=ma if is_loud else None,
        is_bayesian_loud=lambda: is_loud,
        get_model_average_summary_for_model=lambda _model: None,
    )

    write_bayesian_table(report, session)

    return [cell.text for cell in report.document.tables[0].rows[0].cells]


def test_write_bayesian_table_omits_log_posterior_probability_for_loud():
    headers = _bayesian_summary_headers(is_loud=True)

    assert "Unnormalized Log Posterior Probability" not in headers
    assert headers == [
        "Model",
        "Prior Weights",
        "Posterior Weights",
        "BMDL",
        "BMD",
        "BMDU",
        "Scaled Residual at Control",
        "Scaled Residual near BMD",
    ]


def test_write_bayesian_table_keeps_log_posterior_probability_for_non_loud():
    headers = _bayesian_summary_headers(is_loud=False)

    assert "Unnormalized Log Posterior Probability" in headers


def test_ma_model_bmd_triplet_guard_conditions():
    model = SimpleNamespace(settings=SimpleNamespace(alpha=0.1))
    assert _ma_model_bmd_triplet(SimpleNamespace(model_average=None), model) is None

    ma = SimpleNamespace(has_results=False)
    assert _ma_model_bmd_triplet(SimpleNamespace(model_average=ma), model) is None

    ma = SimpleNamespace(has_results=True, models=[], results=SimpleNamespace(model_bmd_dist=[]))
    assert _ma_model_bmd_triplet(SimpleNamespace(model_average=ma), model) is None


def test_ma_model_bmd_triplet_ignores_invalid_draws_and_calculates_quantiles():
    model = SimpleNamespace(settings=SimpleNamespace(alpha=0.25))
    ma = SimpleNamespace(
        has_results=True,
        models=[model],
        results=SimpleNamespace(model_bmd_dist=[[float("nan"), float("inf")]]),
    )
    session = SimpleNamespace(model_average=ma)
    assert _ma_model_bmd_triplet(session, model) is None

    ma.results.model_bmd_dist = [[1.0, 2.0, 3.0, 4.0]]
    assert _ma_model_bmd_triplet(session, model) == (1.75, 2.5, 3.25)
