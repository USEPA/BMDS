import numpy as np
import pandas as pd
import pytest

from pybmds.datasets.transforms.polyk import PolyKAdjustment


@pytest.fixture
def polyk_data(data_path) -> pd.DataFrame:
    return pd.read_csv(data_path / "datasets/transforms/polyk.csv")


class TestAdjustment:
    def test_calculations(self, polyk_data):
        analysis = PolyKAdjustment(
            doses=polyk_data.dose.tolist(),
            day=polyk_data.day.tolist(),
            has_tumor=polyk_data.has_tumor.tolist(),
        )

        adj_df = analysis.adjusted_data
        assert adj_df.day.max() == 734
        assert adj_df.shape == (200, 4)
        assert adj_df.columns.values.tolist() == ["dose", "day", "has_tumor", "adj_n"]

        sum_df = analysis.summary
        assert sum_df.shape == (4, 6)
        assert sum_df.columns.values.tolist() == [
            "dose",
            "n",
            "adj_n",
            "incidence",
            "proportion",
            "adj_proportion",
        ]
        assert np.allclose(sum_df.adj_proportion, [0.1414, 0.2836, 0.5700, 0.6643], atol=1e-4)

    def test_calc_duration_change(self, polyk_data):
        # fixing max_day to a 2 yr cancer bioassay duration
        analysis = PolyKAdjustment(
            doses=polyk_data.dose.tolist(),
            day=polyk_data.day.tolist(),
            has_tumor=polyk_data.has_tumor.tolist(),
            max_day=730,
        )
        assert np.allclose(
            analysis.summary.adj_proportion, [0.1404, 0.2821, 0.5673, 0.6616], atol=1e-4
        )

    def test_reporting(self, polyk_data, data_path, rewrite_data_files):
        analysis = PolyKAdjustment(
            doses=polyk_data.dose.tolist(),
            day=polyk_data.day.tolist(),
            has_tumor=polyk_data.has_tumor.tolist(),
        )

        # excel
        xlsx = analysis.to_excel()

        # docx
        docx = analysis.to_docx()

        if rewrite_data_files:
            (data_path / "reports/data-transform-polyk.xlsx").write_bytes(xlsx.getvalue())
            docx.save(data_path / "reports/data-transform-polyk.docx")
