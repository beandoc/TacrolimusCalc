from pydantic import BaseModel
from typing import List, Optional
from datetime import date, datetime

class ClinicalEventBase(BaseModel):
    datetime: datetime
    dose: Optional[float] = None
    level: Optional[float] = None

class ClinicalEventCreate(ClinicalEventBase):
    pass

class ClinicalEvent(ClinicalEventBase):
    id: int
    patient_mrn: str

    class Config:
        from_attributes = True

class PatientBase(BaseModel):
    mrn: str
    transplant_date: date
    weight: float
    hematocrit: float = 35.0
    mpa: int = 1
    bioassay: int = 1
    genotype: str = "unknown"
    transplant_type: Optional[str] = None
    albumin: float = 3.5
    bilirubin: float = 1.0
    inhibitor: str = "none"

class PatientCreate(PatientBase):
    pass

class Patient(PatientBase):
    events: List[ClinicalEvent] = []

    class Config:
        from_attributes = True
