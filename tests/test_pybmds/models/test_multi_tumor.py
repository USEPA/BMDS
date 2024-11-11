import numpy as np

import pybmds


class TestMultitumor:
    def test_execute(self, mt_datasets, rewrite_data_files, data_path):
        # Check base case where all datasets can be modeled
        degrees = [3, 0, 0]
        session = pybmds.Multitumor(
            mt_datasets, degrees=degrees, id=1, name="test", description="hello"
        )
        session.execute()

        # check text report
        text = session.text()
        assert len(text) > 0

        # check serialization
        session2 = session.serialize().deserialize()
        assert session.to_dict() == session2.to_dict()

        # check that individual models has slope factor and shown in text output
        assert session.models[0][0].results.slope_factor > 0
        assert "Slope Factor" in session.models[0][0].text()

        # dataframe
        df = session.to_df()
        assert "slope_factor" in df.columns
        df = session.params_df()
        df = session.datasets_df()

        # docx
        docx = session.to_docx(all_models=True, bmd_cdf_table=True)

        if rewrite_data_files:
            (data_path / "reports/multitumor.txt").write_text(text)
            df.to_excel(data_path / "reports/multitumor.xlsx", index=False)
            docx.save(data_path / "reports/multitumor.docx")

    def test_one_no_recommend(self, rewrite_data_files, data_path):
        # Check case where one dataset can be modeled, and one has no model recommendation
        datasets = [
            pybmds.DichotomousDataset(
                doses=[0, 25, 75, 125, 200],
                ns=[20, 20, 20, 20, 20],
                incidences=[0, 0, 1, 7, 11],
            ),
            pybmds.DichotomousDataset(
                doses=[0, 50, 100, 200, 400],
                ns=[100, 100, 100, 100, 100],
                incidences=[1, 68, 78, 88, 98],
            ),
        ]
        session = pybmds.Multitumor(datasets)
        session.execute()
        assert session.results.selected_model_indexes == [0, None]
        assert np.isclose(session.results.bmd, 37.69, atol=0.01)
        assert session.results.valid_result is True

        # check we can build reports w/ none
        text = session.text()
        df = session.to_df()
        docx = session.to_docx()

        if rewrite_data_files:
            (data_path / "reports/multitumor-no-selection.txt").write_text(text)
            df.to_excel(data_path / "reports/multitumor-no-selection.xlsx", index=False)
            docx.save(data_path / "reports/multitumor-no-selection.docx")

    def test_all_no_recommend(self):
        # Check case all datasets have no model recommendation
        datasets = [
            pybmds.DichotomousDataset(
                doses=[0, 50, 100, 200, 400],
                ns=[100, 100, 100, 100, 100],
                incidences=[1, 68, 78, 88, 98],
            )
        ]
        session = pybmds.Multitumor(datasets)
        session.execute()

        assert session.results.selected_model_indexes == [None]
        assert session.results.valid_result is False
