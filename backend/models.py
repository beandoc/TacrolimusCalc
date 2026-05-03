from sqlalchemy import Column, Integer, String, Float, DateTime, ForeignKey, Date
from sqlalchemy.orm import relationship
from .database import Base

class Patient(Base):
    __tablename__ = "patients"

    mrn = Column(String, primary_key=True, index=True)
    transplant_date = Column(Date, nullable=False)
    weight = Column(Float, nullable=False)
    hematocrit = Column(Float, default=35.0)
    mpa = Column(Integer, default=1)
    bioassay = Column(Integer, default=1)
    genotype = Column(String, default="unknown")
    transplant_type = Column(String)
    albumin = Column(Float, default=3.5)
    bilirubin = Column(Float, default=1.0)
    inhibitor = Column(String, default="none")
    
    # ML Foundation: Immunological Risk Baseline
    hla_mismatch = Column(Integer, nullable=True) # 0-6
    baseline_pra = Column(Float, nullable=True)   # % PRA
    
    # ML Foundation: Advanced SOTA Predictors (XGBoost 30-Day Rejection Model)
    recipient_age = Column(Integer, nullable=True)
    recipient_bsa = Column(Float, nullable=True)  # Dubois formula
    retransplant = Column(Integer, default=0)     # 0=False, 1=True
    donor_age = Column(Integer, nullable=True)
    donor_bsa = Column(Float, nullable=True)
    donor_vasoactive_drugs = Column(Integer, default=0) # 0=False, 1=True
    induction_therapy = Column(String, nullable=True)   # e.g., Thymoglobulin, Basiliximab
    dgf = Column(Integer, default=0)              # Delayed Graft Function (0=No, 1=Yes)
    cold_ischemia_time = Column(Float, nullable=True)   # Hours

    events = relationship("ClinicalEvent", back_populates="patient", cascade="all, delete-orphan")
    outcomes = relationship("ClinicalOutcome", back_populates="patient", cascade="all, delete-orphan")

class ClinicalEvent(Base):
    __tablename__ = "clinical_events"

    id = Column(Integer, primary_key=True, index=True)
    patient_mrn = Column(String, ForeignKey("patients.mrn"))
    datetime = Column(DateTime, nullable=False)
    dose = Column(Float, nullable=True)
    level = Column(Float, nullable=True)
    
    # ML Foundation: Biomarkers
    creatinine = Column(Float, nullable=True)
    wbc = Column(Float, nullable=True)
    crp = Column(Float, nullable=True)

    patient = relationship("Patient", back_populates="events")

class ClinicalOutcome(Base):
    __tablename__ = "clinical_outcomes"

    id = Column(Integer, primary_key=True, index=True)
    patient_mrn = Column(String, ForeignKey("patients.mrn"))
    diagnosis_date = Column(DateTime, nullable=False)
    diagnosis_type = Column(String, nullable=False) # e.g. ACR, AMR, CNI Toxicity, Infection, BKVAN
    biopsy_proven = Column(Integer, default=0) # 0=Clinical, 1=Biopsy-proven
    notes = Column(String, nullable=True)

    patient = relationship("Patient", back_populates="outcomes")
