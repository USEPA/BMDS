from types import SimpleNamespace

from pybmds.reporting.styling import Report, write_bayesian_table, write_dataset_table


def test_write_dataset_table(cdataset, cidataset, ddataset, nd_dataset):
    report = Report.build_default()
    write_dataset_table(report, cdataset, True)
    write_dataset_table(report, cdataset, False)
    write_dataset_table(report, cidataset, True)
    write_dataset_table(report, cidataset, False)
    write_dataset_table(report, ddataset, True)
    write_dataset_table(report, ddataset, False)
    write_dataset_table(report, nd_dataset, False)


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
