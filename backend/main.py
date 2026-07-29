import sys
import os
from fastapi import FastAPI, Depends, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy import inspect, text
from sqlalchemy.orm import Session
from typing import List

try:
    from . import models, schemas
    from .database import engine, SessionLocal
except ImportError:
    sys.path.insert(0, os.path.dirname(__file__))
    import models
    import schemas
    from database import engine, SessionLocal

# Create database tables
models.Base.metadata.create_all(bind=engine)


def _sync_missing_columns():
    # create_all() only creates tables that don't exist yet — it never ALTERs
    # a table that's already there. So a column added to a model after the
    # table was first created (e.g. Patient.height) silently never reaches a
    # live Postgres database, and every query touching that table then 500s
    # with psycopg2.errors.UndefinedColumn. This adds any such missing
    # columns on startup so model changes stay in sync with deployed DBs
    # without needing a separate migration tool.
    inspector = inspect(engine)
    for table in models.Base.metadata.sorted_tables:
        if not inspector.has_table(table.name):
            continue
        existing_cols = {c["name"] for c in inspector.get_columns(table.name)}
        for col in table.columns:
            if col.name in existing_cols:
                continue
            col_type = col.type.compile(dialect=engine.dialect)
            try:
                with engine.begin() as conn:
                    conn.execute(text(f'ALTER TABLE {table.name} ADD COLUMN "{col.name}" {col_type}'))
            except Exception as e:
                print(f"schema sync: failed to add {table.name}.{col.name}: {e}")


_sync_missing_columns()

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
        # Update existing with only the fields that were actually set
        update_data = patient.model_dump(exclude_unset=True)
        for var, value in update_data.items():
            setattr(db_patient, var, value)
    else:
        # Create new
        db_patient = models.Patient(**patient.model_dump())
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
        db_event = models.ClinicalEvent(**event.model_dump(), patient_mrn=mrn)
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
        db_outcome = models.ClinicalOutcome(**outcome.model_dump(), patient_mrn=mrn)
        db.add(db_outcome)
        new_outcomes.append(db_outcome)
        
    db.commit()
    return db.query(models.ClinicalOutcome).filter(models.ClinicalOutcome.patient_mrn == mrn).all()

@app.get("/api/patients/{mrn}/outcomes", response_model=List[schemas.ClinicalOutcome])
def read_outcomes(mrn: str, db: Session = Depends(get_db)):
    outcomes = db.query(models.ClinicalOutcome).filter(models.ClinicalOutcome.patient_mrn == mrn).all()
    return outcomes
