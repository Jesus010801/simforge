"""Source assertions about chemistry; no equivalence or chemistry execution."""
from typing import Annotated, Literal
from pydantic import Field, model_validator
from .entities import Entity
from .identity import _ChemicalRef, _ImportProvenance
from .serialization import _DomainModel, _Items, _Text, _Symbol


class _ChemicalRepresentation(_DomainModel):
    representation_type: _Symbol
    value: _Text
    specificity: Literal['PARTIAL', 'SPECIFIED']
    provenance: _ImportProvenance


class ChemicalIdentity(Entity):
    entity_type: Literal['chemical_identity'] = 'chemical_identity'
    label: _Text | None = None
    representations: Annotated[_Items[_ChemicalRepresentation], Field(min_length=1)]
    stereochemistry_status: Literal['SPECIFIED', 'UNSPECIFIED', 'NOT_APPLICABLE']
    isotope_status: Literal['SPECIFIED', 'UNSPECIFIED', 'NOT_APPLICABLE']
    _unordered_fields = ('representations',)


class ChemicalSpecies(Entity):
    entity_type: Literal['chemical_species'] = 'chemical_species'
    chemical_identity: _ChemicalRef
    representations: Annotated[_Items[_ChemicalRepresentation], Field(min_length=1)]
    formal_charge: int | None
    protonation_state: _Text | None
    tautomer_state: _Text | None
    description_status: Literal['PARTIAL', 'SPECIFIED']
    _unordered_fields = ('representations',)

    @model_validator(mode='after')
    def _specificity(self):
        if self.sf_id == self.chemical_identity.entity_id:
            raise ValueError('species and parent identity must have distinct IDs')
        if self.description_status == 'SPECIFIED' and any(v is None for v in (self.formal_charge, self.protonation_state, self.tautomer_state)):
            raise ValueError('specified description requires explicit species details')
        return self
