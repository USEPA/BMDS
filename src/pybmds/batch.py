import json
import os
import zipfile
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from io import BytesIO
from pathlib import Path
from typing import NamedTuple, Self

import pandas as pd
from tqdm import tqdm

from .datasets.base import DatasetType
from .models.multi_tumor import Multitumor
from .reporting.styling import Report, write_citation
from .session import Session


class BatchResponse(NamedTuple):
    success: bool
    content: dict | list[dict]


class BatchBase:
    pass


class BatchSession(BatchBase):
    def __init__(self, sessions: list[Session] | None = None):
        if sessions is None:
            sessions = []
        self.sessions: list[Session] = sessions
        self.errors = []

    def df_summary(self) -> pd.DataFrame:
        dfs = [
            session.to_df(
                extras=dict(
                    session_index=idx,
                    session_id=session.id,
                    session_name=session.name,
                    session_description=session.description,
                ),
                clean=False,
            )
            for idx, session in enumerate(self.sessions)
        ]
        return pd.concat(dfs).dropna(axis=1, how="all").fillna("")

    def df_dataset(self) -> pd.DataFrame:
        data: list[dict] = []
        for idx, session in enumerate(self.sessions):
            data.extend(
                session.dataset.rows(
                    extras=dict(
                        session_index=idx,
                        session_id=session.id,
                        session_name=session.name,
                        session_description=session.description,
                    )
                )
            )
        return pd.DataFrame(data=data)

    def df_params(self) -> pd.DataFrame:
        data: list[dict] = []
        for idx, session in enumerate(self.sessions):
            for model_index, model in enumerate(session.models):
                if model.has_results:
                    data.extend(
                        model.results.parameters.rows(
                            extras=dict(
                                session_index=idx,
                                session_id=session.id,
                                session_name=session.name,
                                session_description=session.description,
                                model_index=model_index,
                                model_name=model.name(),
                            )
                        )
                    )
        return pd.DataFrame(data=data)

    def to_excel(self, path: Path | None = None) -> Path | BytesIO:
        f: Path | BytesIO = path or BytesIO()
        with pd.ExcelWriter(f) as writer:
            data = {
                "summary": self.df_summary(),
                "datasets": self.df_dataset(),
                "parameters": self.df_params(),
            }
            for name, df in data.items():
                df.to_excel(writer, sheet_name=name, index=False)
        return f

    def to_docx(
        self,
        report: Report | None = None,
        header_level: int = 1,
        citation: bool = True,
        dataset_format_long: bool = True,
        all_models: bool = False,
        bmd_cdf_table: bool = False,
        session_inputs_table: bool = False,
    ):
        """Append each session to a single document

        Args:
            report (Report, optional): A Report object, or None to use default.
            header_level (int, optional): Starting header level. Defaults to 1.
            citation (bool, default True): Include citation
            dataset_format_long (bool, default True): long or wide dataset table format
            all_models (bool, default False):  Show all models, not just selected
            bmd_cdf_table (bool, default False): Export BMD CDF table
            session_inputs_table (bool, default False): Write an inputs table for a session,
                assuming a single model's input settings are representative of all models in a
                session, which may not always be true

        Returns:
            A python docx.Document object with content added.
        """
        if report is None:
            report = Report.build_default()

        for session in self.sessions:
            session.to_docx(
                report,
                header_level=header_level,
                citation=False,
                dataset_format_long=dataset_format_long,
                all_models=all_models,
                bmd_cdf_table=bmd_cdf_table,
                session_inputs_table=session_inputs_table,
            )

        if citation and len(self.sessions) > 0:
            write_citation(report, header_level=header_level)

        return report.document

    def serialize(self) -> str:
        """Export Session into a json format which can be saved and loaded.

        Returns:
            str: A JSON string
        """
        return json.dumps([session.to_dict() for session in self.sessions])

    @classmethod
    def execute(
        cls,
        datasets: list[DatasetType],
        runner: Callable[[DatasetType], BatchResponse],
        nprocs: int | None = None,
    ) -> Self:
        """Execute sessions using multiple processors.

        Args:
            datasets (list[DatasetType]): The datasets to execute
            runner (Callable[dataset] -> BatchResponse): The method which executes a session
            nprocs (Optional[int]): the number of processors to use; defaults to N-1. If 1 is
                specified; the batch session is called linearly without a process pool

        Returns:
            A BatchSession with sessions executed.
        """
        if nprocs is None:
            nprocs = max(os.cpu_count() - 1, 1)

        results: list[BatchResponse] = []
        if nprocs == 1:
            # run without a ProcessPoolExecutor; useful for debugging
            for dataset in tqdm(datasets, desc="Executing..."):
                results.append(runner(dataset))
        else:
            # adapted from https://gist.github.com/alexeygrigorev/79c97c1e9dd854562df9bbeea76fc5de
            with ProcessPoolExecutor(max_workers=nprocs) as executor:
                with tqdm(total=len(datasets), desc="Executing...") as progress:
                    futures = []
                    for dataset in datasets:
                        future = executor.submit(runner, dataset)
                        future.add_done_callback(lambda p: progress.update())
                        futures.append(future)

                    for future in futures:
                        results.append(future.result())

        batch = cls()
        for result in tqdm(results, desc="Building batch..."):
            if result.success:
                if isinstance(result.content, list):
                    for item in result.content:
                        batch.sessions.append(Session.from_serialized(item))
                else:
                    batch.sessions.append(Session.from_serialized(result.content))
            else:
                batch.errors.append(result.content)

        return batch

    @classmethod
    def deserialize(cls, data: str) -> Self:
        """Load serialized batch export into a batch session.

        Args:
            data (str): A JSON export generated from the `BatchSession.serialize` method.
        """
        sessions_data = json.loads(data)
        sessions = [Session.from_serialized(session) for session in sessions_data]
        return cls(sessions=sessions)

    @classmethod
    def load(cls, archive: Path) -> Self:
        """Load a Session from a compressed zipfile

        Args:
            fn (Path): The zipfile path

        Returns:
            BatchSession: An instance of this session
        """
        with zipfile.ZipFile(archive) as zf:
            with zf.open("data.json") as f:
                data = f.read()
        return BatchSession.deserialize(data)

    def save(self, archive: Path):
        """Save Session to a compressed zipfile

        Args:
            fn (Path): The zipfile path
        """
        with zipfile.ZipFile(
            archive, mode="w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
        ) as zf:
            zf.writestr("data.json", data=self.serialize())


class MultitumorBatch(BatchBase):
    def __init__(self, sessions: list[Multitumor] | None = None):
        if sessions is None:
            sessions = []
        self.sessions: list[Multitumor] = sessions
        self.errors = []

    def to_docx(
        self,
        report: Report | None = None,
        header_level: int = 1,
        citation: bool = True,
        **kw,
    ):
        """Append each session to a single document

        Args:
            report (Report, optional): A Report object, or None to use default.
            header_level (int, optional): Starting header level. Defaults to 1.
            citation (bool, default True): Include citation

        Returns:
            A python docx.Document object with content added.
        """
        if report is None:
            report = Report.build_default()

        for session in self.sessions:
            session.to_docx(
                report,
                header_level=header_level,
                citation=False,
            )

        if citation and len(self.sessions) > 0:
            write_citation(report, header_level=header_level)

        return report.document

    def serialize(self) -> str:
        return json.dumps([session.to_dict() for session in self.sessions])

    @classmethod
    def deserialize(cls, data: str) -> Self:
        sessions_data = json.loads(data)
        sessions = [Multitumor.from_serialized(session) for session in sessions_data]
        return cls(sessions=sessions)

    def df_summary(self) -> pd.DataFrame:
        dfs = [
            session.to_df(
                extras=dict(session_index=idx),
                clean=False,
            )
            for idx, session in enumerate(self.sessions)
        ]
        return pd.concat(dfs).dropna(axis=1, how="all").fillna("")

    def df_dataset(self) -> pd.DataFrame:
        dfs = [
            session.datasets_df(
                extras=dict(session_index=idx),
            )
            for idx, session in enumerate(self.sessions)
        ]
        return pd.concat(dfs).dropna(axis=1, how="all").fillna("")

    def df_params(self) -> pd.DataFrame:
        dfs = [
            session.params_df(
                extras=dict(session_index=idx),
            )
            for idx, session in enumerate(self.sessions)
        ]
        return pd.concat(dfs).dropna(axis=1, how="all").fillna("")

    def to_excel(self, path: Path | None = None) -> Path | BytesIO:
        f: Path | BytesIO = path or BytesIO()
        with pd.ExcelWriter(f) as writer:
            data = {
                "summary": self.df_summary(),
                "datasets": self.df_dataset(),
                "parameters": self.df_params(),
            }
            for name, df in data.items():
                df.to_excel(writer, sheet_name=name, index=False)
        return f
