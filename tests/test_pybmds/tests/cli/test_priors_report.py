from pybmds.cli.priors_report import create_report


def test_create_report():
    # ensure that the report method works
    report = create_report().getvalue()
    assert report.startswith("# BMDS priors report")
