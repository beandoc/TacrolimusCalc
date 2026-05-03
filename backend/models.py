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

    events = relationship("ClinicalEvent", back_populates="patient", cascade="all, delete-orphan")

class ClinicalEvent(Base):
    __tablename__ = "clinical_events"

    id = Column(Integer, primary_key=True, index=True)
    patient_mrn = Column(String, ForeignKey("patients.mrn"))
    datetime = Column(DateTime, nullable=False)
    dose = Column(Float, nullable=True)
    level = Column(Float, nullable=True)

    patient = relationship("Patient", back_populates="events")
