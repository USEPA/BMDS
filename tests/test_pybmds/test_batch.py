import tempfile
from pathlib import Path

import pybmds
from pybmds.batch import BatchResponse, BatchSession, MultitumorBatch
from pybmds.session import Session


def _batch_run(ds):
    # added to root of file so that we can pickle for multiprocessing
    sess = Session(dataset=ds)
    sess.add_model(pybmds.Models.Logistic)
    sess.execute()
    return BatchResponse(success=True, content=[sess.to_dict()])


class TestBatchSession:
    def test_execute(self, ddataset2):
        batch = BatchSession.execute([ddataset2], _batch_run, nprocs=1)
        assert len(batch.session) == 1

        batch = BatchSession.execute([ddataset2, ddataset2], _batch_run, nprocs=2)
        assert len(batch.session) == 2

    def test_exports_dichotomous(self, ddataset2, rewrite_data_files, data_path):
        datasets = [ddataset2]
        batch = BatchSession()
        for dataset in datasets:
            session = Session(dataset=dataset)
            session.add_default_models()
            session.execute_and_recommend()
            batch.session.append(session)

            session = Session(dataset=dataset)
            session.add_default_bayesian_models()
            session.execute()
            batch.session.append(session)

        # check serialization/deserialization
        data = batch.serialize()
        batch2 = batch.deserialize(data)
        assert len(batch2.session) == len(batch.session)

        # check zip
        zf = Path(tempfile.NamedTemporaryFile().name)
        try:
            # save
            batch.save(zf)
            assert zf.exists()
            # unsave
            batch3 = BatchSession.load(zf)
            assert len(batch3.session) == 2
        finally:
            zf.unlink()

        # check exports
        excel = batch.to_excel()
        docx = batch.to_docx()

        if rewrite_data_files:
            (data_path / "reports/batch-dichotomous.xlsx").write_bytes(excel.getvalue())
            docx.save(data_path / "reports/batch-dichotomous.docx")

    def test_exports_continuous(self, cdataset2, cidataset, rewrite_data_files, data_path):
        datasets = [cdataset2, cidataset]
        batch = BatchSession()
        for dataset in datasets:
            session = pybmds.Session(dataset=dataset)
            session.add_model(pybmds.Models.Power)
            session.execute_and_recommend()
            batch.session.append(session)

        # check serialization/deserialization
        data = batch.serialize()
        batch2 = batch.deserialize(data)
        assert len(batch2.session) == len(batch.session)

        # check exports
        excel = batch.to_excel()
        docx = batch.to_docx()

        if rewrite_data_files:
            (data_path / "reports/batch-continuous.xlsx").write_bytes(excel.getvalue())
            docx.save(data_path / "reports/batch-continuous.docx")


class TestMultitumorBatch:
    def test_exports(self, mt_datasets, rewrite_data_files, data_path):
        session = pybmds.Multitumor(
            mt_datasets, degrees=[0] * len(mt_datasets), id=1, name="test", description="hello"
        )
        session.execute()
        batch = MultitumorBatch(sessions=[session])

        # check serialization/deserialization
        data = batch.serialize()
        batch2 = batch.deserialize(data)
        assert len(batch2.session) == len(batch.session)

        # check exports
        excel = batch.to_excel()
        docx = batch.to_docx()

        if rewrite_data_files:
            (data_path / "reports/batch-multitumor.xlsx").write_bytes(excel.getvalue())
            docx.save(data_path / "reports/batch-multitumor.docx")
