from pydantic import BaseModel
from typing import List, Optional
from datetime import date, datetime

class ClinicalEventBase(BaseModel):
    datetime: datetime
    dose: Optional[float] = None
    level: Optional[float] = None
    creatinine: Optional[float] = None
    wbc: Optional[float] = None
    crp: Optional[float] = None

class ClinicalEventCreate(ClinicalEventBase):
    pass

class ClinicalEvent(ClinicalEventBase):
    id: int
    patient_mrn: str

    class Config:
        from_attributes = True

class ClinicalOutcomeBase(BaseModel):
    diagnosis_date: datetime
    diagnosis_type: str
    biopsy_proven: int = 0
    notes: Optional[str] = None

class ClinicalOutcomeCreate(ClinicalOutcomeBase):
    pass

class ClinicalOutcome(ClinicalOutcomeBase):
    id: int
    patient_mrn: str

    class Config:
        from_attributes = True

class PatientBase(BaseModel):
    mrn: str
    transplant_date: date
    weight: float
    height: Optional[float] = None
    hematocrit: float = 35.0
    mpa: int = 1
    bioassay: int = 1
    genotype: str = "unknown"
    transplant_type: Optional[str] = None
    albumin: float = 3.5
    bilirubin: float = 1.0
    inhibitor: str = "none"
    hla_mismatch: Optional[int] = None
    baseline_pra: Optional[float] = None
    
    # Advanced SOTA Predictors
    recipient_age: Optional[int] = None
    recipient_bsa: Optional[float] = None
    retransplant: int = 0
    donor_age: Optional[int] = None
    donor_bsa: Optional[float] = None
    donor_vasoactive_drugs: int = 0
    induction_therapy: Optional[str] = None
    dgf: int = 0
    cold_ischemia_time: Optional[float] = None

class PatientCreate(PatientBase):
    pass

class Patient(PatientBase):
    events: List[ClinicalEvent] = []
    outcomes: List[ClinicalOutcome] = []

    class Config:
        from_attributes = True
