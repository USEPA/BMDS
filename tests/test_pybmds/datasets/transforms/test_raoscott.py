import numpy as np
import pytest

from pybmds.datasets import DichotomousDataset
from pybmds.datasets.transforms.raoscott import RaoScott, Species


@pytest.fixture
def dataset() -> DichotomousDataset:
    return DichotomousDataset(
        doses=[0, 7, 35, 100, 175, 350, 500],
        ns=[470, 211, 232, 220, 241, 237, 166],
        incidences=[11, 6, 2, 7, 14, 39, 57],
    )


class TestRaoScott:
    def test_calculations(self, dataset):
        analysis = RaoScott(dataset=dataset, species=Species.rat)

        assert np.allclose(
            analysis.df.n_adjusted,
            [284.16, 119.16, 198.92, 119.45, 105.98, 72.36, 39.15],
            atol=0.01,
        )
        assert np.allclose(
            analysis.df.incidence_adjusted,
            [6.65, 3.39, 1.71, 3.8, 6.16, 11.91, 13.44],
            atol=0.01,
        )

        assert np.allclose(
            analysis.df.incidence / analysis.df.n,
            analysis.df.incidence_adjusted / analysis.df.n_adjusted,
        )

    def test_reporting(self, dataset, data_path, rewrite_data_files):
        analysis = RaoScott(dataset=dataset, species=Species.rat)

        xlsx = analysis.to_excel()
        docx = analysis.to_docx()

        if rewrite_data_files:
            (data_path / "reports/data-transform-rao-scott.xlsx").write_bytes(xlsx.getvalue())
            docx.save(data_path / "reports/data-transform-rao-scott.docx")
