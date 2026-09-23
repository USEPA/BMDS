import sys
from io import StringIO

from pybmds.cli import priors_report
from pybmds.cli.priors_report import create_report


def test_create_report():
    # ensure that the report method works
    report = create_report().getvalue()
    assert report.startswith("# BMDS priors report")


def test_main_writes_report_to_stdout(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["priors_report"])
    monkeypatch.setattr(priors_report, "create_report", lambda: StringIO("report text"))

    priors_report.main()

    assert capsys.readouterr().out == "report text"


def test_main_writes_report_to_file(monkeypatch, capsys):
    written = {}
    monkeypatch.setattr(sys, "argv", ["priors_report", "priors.md"])
    monkeypatch.setattr(priors_report, "create_report", lambda: StringIO("report text"))
    monkeypatch.setattr(
        priors_report.Path,
        "write_text",
        lambda self, text: written.update(path=self, text=text),
    )

    priors_report.main()

    assert written["path"].name == "priors.md"
    assert written["text"] == "report text"
    assert "Writing output to:" in capsys.readouterr().out
