from pydantic import BaseModel, Field

from .models.base import BmdModel


class SelectedModelSchema(BaseModel):
    bmds_model_index: int | None = Field(default=None, alias="model_index")
    notes: str = ""

    def deserialize(self, session) -> "SelectedModel":
        return SelectedModel(session, model_index=self.bmds_model_index, notes=self.notes)


class SelectedModel:
    def __init__(self, session, model_index: int | None = None, notes: str = ""):
        self.session = session
        self.model_index = model_index
        self.notes = notes

    def select(self, model: BmdModel | None, notes: str):
        self.model_index = self.session.models.index(model) if model is not None else None
        self.notes = notes

    @property
    def model(self) -> BmdModel | None:
        """Returns the selected model if one exists, else None"""
        if self.model_index is not None:
            return self.session.models[self.model_index]
        return None

    @property
    def no_model_selected(self) -> bool:
        """Check if no model selected was deliberate or just undefined.

        Assumes that if a user has provided notes then that is a clear decision that not model was
        selected; otherwise it's likely that a selection is just undefined.
        """
        return self.model_index is None and isinstance(self.notes, str)

    def serialize(self) -> SelectedModelSchema:
        return SelectedModelSchema(model_index=self.model_index, notes=self.notes)

    def update_record(self, d: dict, index: int) -> None:
        """Update data record for a tabular-friendly export"""
        is_selected = self.model_index == index
        d.update(
            selected=is_selected,
            selected_notes=self.notes if is_selected else None,
        )
