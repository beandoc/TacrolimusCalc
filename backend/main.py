from fastapi import FastAPI, Depends, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy.orm import Session
from typing import List

from . import models, schemas
from .database import engine, SessionLocal

# Create database tables
models.Base.metadata.create_all(bind=engine)

app = FastAPI(title="Tacrolimus API")

# Configure CORS for Vercel integration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], # In production, restrict to your specific Vercel URL
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Dependency to get DB session
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

@app.get("/api/patients", response_model=List[schemas.Patient])
def read_patients(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    patients = db.query(models.Patient).offset(skip).limit(limit).all()
    return patients

@app.post("/api/patients", response_model=schemas.Patient)
def create_patient(patient: schemas.PatientCreate, db: Session = Depends(get_db)):
    db_patient = db.query(models.Patient).filter(models.Patient.mrn == patient.mrn).first()
    if db_patient:
        # Update existing
        for var, value in vars(patient).items():
            setattr(db_patient, var, value)
    else:
        # Create new
        db_patient = models.Patient(**patient.dict())
        db.add(db_patient)
    db.commit()
    db.refresh(db_patient)
    return db_patient

@app.get("/api/patients/{mrn}", response_model=schemas.Patient)
def read_patient(mrn: str, db: Session = Depends(get_db)):
    db_patient = db.query(models.Patient).filter(models.Patient.mrn == mrn).first()
    if db_patient is None:
        raise HTTPException(status_code=404, detail="Patient not found")
    return db_patient

@app.post("/api/patients/{mrn}/events", response_model=List[schemas.ClinicalEvent])
def create_events(mrn: str, events: List[schemas.ClinicalEventCreate], db: Session = Depends(get_db)):
    db_patient = db.query(models.Patient).filter(models.Patient.mrn == mrn).first()
    if not db_patient:
        raise HTTPException(status_code=404, detail="Patient not found")
    
    # Optional: Clear old events to sync completely, or just append
    db.query(models.ClinicalEvent).filter(models.ClinicalEvent.patient_mrn == mrn).delete()
    
    new_events = []
    for event in events:
        db_event = models.ClinicalEvent(**event.dict(), patient_mrn=mrn)
        db.add(db_event)
        new_events.append(db_event)
        
    db.commit()
    
    # Return updated list
    return db.query(models.ClinicalEvent).filter(models.ClinicalEvent.patient_mrn == mrn).all()

@app.get("/api/patients/{mrn}/events", response_model=List[schemas.ClinicalEvent])
def read_events(mrn: str, db: Session = Depends(get_db)):
    events = db.query(models.ClinicalEvent).filter(models.ClinicalEvent.patient_mrn == mrn).all()
    return events

@app.post("/api/patients/{mrn}/outcomes", response_model=List[schemas.ClinicalOutcome])
def create_outcomes(mrn: str, outcomes: List[schemas.ClinicalOutcomeCreate], db: Session = Depends(get_db)):
    db_patient = db.query(models.Patient).filter(models.Patient.mrn == mrn).first()
    if not db_patient:
        raise HTTPException(status_code=404, detail="Patient not found")
    
    db.query(models.ClinicalOutcome).filter(models.ClinicalOutcome.patient_mrn == mrn).delete()
    
    new_outcomes = []
    for outcome in outcomes:
        db_outcome = models.ClinicalOutcome(**outcome.dict(), patient_mrn=mrn)
        db.add(db_outcome)
        new_outcomes.append(db_outcome)
        
    db.commit()
    return db.query(models.ClinicalOutcome).filter(models.ClinicalOutcome.patient_mrn == mrn).all()

@app.get("/api/patients/{mrn}/outcomes", response_model=List[schemas.ClinicalOutcome])
def read_outcomes(mrn: str, db: Session = Depends(get_db)):
    outcomes = db.query(models.ClinicalOutcome).filter(models.ClinicalOutcome.patient_mrn == mrn).all()
    return outcomes
