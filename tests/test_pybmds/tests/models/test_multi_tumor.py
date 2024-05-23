from pybmds.models.multi_tumor import Multitumor


class TestMultitumor:
    def test_execute(self, mt_datasets, rewrite_data_files, data_path):
        degrees = [3, 0, 0]
        session = Multitumor(mt_datasets, degrees=degrees, id=1, name="test", description="hello")
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
        docx = session.to_docx()

        if rewrite_data_files:
            (data_path / "reports/multitumor.txt").write_text(text)
            df.to_excel(data_path / "reports/multitumor.xlsx", index=False)
            docx.save(data_path / "reports/multitumor.docx")
