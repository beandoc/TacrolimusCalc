from pydantic import BaseModel, Field
from typing import List, Optional
from datetime import date, datetime

class ClinicalEventBase(BaseModel):
    datetime: datetime
    dose: Optional[float] = None
    level: Optional[float] = None
    weight: Optional[float] = None
    hematocrit: Optional[float] = None
    creatinine: Optional[float] = None
    wbc: Optional[float] = None
    crp: Optional[float] = None
    inhibitor: Optional[str] = None

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
    # 0 = on tacrolimus, 1 = switched to another agent. Optional because rows
    # created before this column exists deserialise as None; consumers must read
    # it as falsy-means-active rather than assuming 0.
    tac_discontinued: Optional[int] = 0
    discontinued_reason: Optional[str] = None
    hla_mismatch: Optional[int] = None
    baseline_pra: Optional[float] = None
    
    # Advanced SOTA Predictors
    recipient_age: Optional[int] = None
    recipient_bsa: Optional[float] = None
    retransplant: Optional[int] = 0
    donor_age: Optional[int] = None
    donor_bsa: Optional[float] = None
    donor_vasoactive_drugs: Optional[int] = 0
    induction_therapy: Optional[str] = None
    dgf: Optional[int] = 0
    cold_ischemia_time: Optional[float] = None

class PatientCreate(PatientBase):
    pass

class PatientSummary(PatientBase):
    """Patient WITHOUT the event/outcome collections.

    The list endpoint fills a name picker and the Load Patient directory, which
    need mrn / transplant_date / weight and nothing else. Returning the full
    `Patient` there serialised every dose and level in the database on each call
    (plus one lazy-load query per patient) — so opening the picker downloaded the
    entire clinical record set.
    """

    class Config:
        from_attributes = True

class Patient(PatientBase):
    events: List[ClinicalEvent] = []
    outcomes: List[ClinicalOutcome] = []

    class Config:
        from_attributes = True

class CenterConfigUpdate(BaseModel):
    # Bounded deliberately. This single number multiplies TVCL for EVERY patient
    # on EVERY device that loads the app, so an out-of-range value is not a
    # cosmetic bug: cl_scalar = 0 gives CL = 0 (no elimination, unbounded
    # predicted concentrations), and a negative value gives negative clearance.
    #
    # The plausible range for a centre-level correction to the Størset prior is
    # narrow — the shipped default is 1.18, and a centre needing more than ±60%
    # has a data problem (assay calibration, dose-log quality) that recalibrating
    # the prior would only paper over.
    cl_scalar: float = Field(gt=0.5, lt=3.0, description="Centre-wide multiplier on TVCL")

class CenterConfig(CenterConfigUpdate):
    updated_at: Optional[datetime] = None

    class Config:
        from_attributes = True
