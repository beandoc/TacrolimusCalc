import os
import secrets
import sys
from datetime import datetime, timezone
from fastapi import Depends, FastAPI, Header, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy import inspect, text
from sqlalchemy.orm import Session
from typing import List, Optional

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

# ---------------------------------------------------------------------------
# CORS
#
# `allow_origins=["*"]` together with `allow_credentials=True` is rejected by
# browsers (the CORS spec forbids the wildcard on credentialed requests), so the
# previous config was not doing what it looked like — while still advertising
# every endpoint to every origin. Set TACRO_ALLOWED_ORIGINS to a comma-separated
# list of the origins that actually serve the app.
#
# The default stays permissive ONLY for local development, where the app is
# routinely opened over file:// (origin "null") against localhost:8000.
# ---------------------------------------------------------------------------
_origins_env = os.getenv("TACRO_ALLOWED_ORIGINS", "").strip()
if _origins_env:
    _allowed_origins = [o.strip() for o in _origins_env.split(",") if o.strip()]
else:
    _allowed_origins = ["*"]
    print(
        "WARNING: TACRO_ALLOWED_ORIGINS is unset — allowing all origins. "
        "Set it to your deployed origin(s) before exposing this API."
    )

app.add_middleware(
    CORSMiddleware,
    allow_origins=_allowed_origins,
    # Auth travels in an explicit header, not a cookie, so credentialed requests
    # are never needed — and leaving this False is what makes a wildcard legal.
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ---------------------------------------------------------------------------
# AUTHENTICATION
#
# This API stores patient-identifiable data: names/MRNs, transplant dates and
# full dose/level histories. It previously had no authentication of any kind, so
# anything that could reach the host could enumerate and overwrite the entire
# clinical database.
#
# A single shared token is the smallest thing that actually closes that hole and
# fits how the tool is used (one transplant unit, a handful of trusted devices).
# It is deliberately NOT a user-identity system — there is no per-user audit
# trail, and it should not be mistaken for one.
#
# Set TACRO_API_TOKEN on the server; the browser sends it as `X-Tacro-Token`,
# stored alongside the backend URL under Settings. If the variable is unset the
# API stays open and says so loudly at startup, so existing local setups keep
# working — but a deployment without it is not protected.
# ---------------------------------------------------------------------------
API_TOKEN = os.getenv("TACRO_API_TOKEN", "").strip()
if not API_TOKEN:
    print(
        "WARNING: TACRO_API_TOKEN is unset — the API is UNAUTHENTICATED and will "
        "serve patient data to any client that can reach it. Set it before deploying."
    )


def require_token(x_tacro_token: Optional[str] = Header(default=None)):
    """Reject requests without the shared token, when one is configured."""
    if not API_TOKEN:
        return
    # compare_digest avoids leaking the token through response-timing.
    if not x_tacro_token or not secrets.compare_digest(x_tacro_token, API_TOKEN):
        raise HTTPException(status_code=401, detail="Invalid or missing API token")


# Dependency to get DB session
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

@app.get("/api/patients", response_model=List[schemas.PatientSummary])
def read_patients(skip: int = 0, limit: int = 500, db: Session = Depends(get_db), _=Depends(require_token)):
    # PatientSummary, not Patient: the full schema carries `events` and
    # `outcomes`, so listing patients serialised every dose and level in the
    # database (plus an N+1 lazy load per patient) just to fill a name picker.
    patients = db.query(models.Patient).offset(skip).limit(limit).all()
    return patients

@app.post("/api/patients", response_model=schemas.Patient)
def create_patient(patient: schemas.PatientCreate, db: Session = Depends(get_db), _=Depends(require_token)):
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
def read_patient(mrn: str, db: Session = Depends(get_db), _=Depends(require_token)):
    db_patient = db.query(models.Patient).filter(models.Patient.mrn == mrn).first()
    if db_patient is None:
        raise HTTPException(status_code=404, detail="Patient not found")
    return db_patient

@app.post("/api/patients/{mrn}/events", response_model=List[schemas.ClinicalEvent])
def create_events(
    mrn: str,
    events: List[schemas.ClinicalEventCreate],
    db: Session = Depends(get_db),
    allow_shrink: bool = Query(
        False,
        description="Permit a write that reduces the number of stored events. "
                    "Required for deletions; blocks accidental history loss.",
    ),
    _=Depends(require_token),
):
    db_patient = db.query(models.Patient).filter(models.Patient.mrn == mrn).first()
    if not db_patient:
        raise HTTPException(status_code=404, detail="Patient not found")

    # This endpoint REPLACES the whole event history, which makes it a data-loss
    # weapon in one specific, entirely realistic sequence: the backend blips
    # while loadPatientData() runs, the client falls back to a stale or empty
    # localStorage copy, the user adds one record, and the autosave that follows
    # overwrites months of history with that fallback. Every add/edit/delete in
    # the UI autosaves, so it needs no unusual action to trigger.
    #
    # Refuse any write that would REDUCE the stored count unless the caller says
    # explicitly that shrinking is intended. Growing or same-size writes (the
    # overwhelming majority) are unaffected.
    existing_count = (
        db.query(models.ClinicalEvent).filter(models.ClinicalEvent.patient_mrn == mrn).count()
    )
    if len(events) < existing_count and not allow_shrink:
        raise HTTPException(
            status_code=409,
            detail=(
                f"Refusing to replace {existing_count} stored event(s) with {len(events)}. "
                "This usually means the client is working from a stale or partial copy. "
                "Reload the patient, or resend with ?allow_shrink=true if the removal is intended."
            ),
        )

    db.query(models.ClinicalEvent).filter(models.ClinicalEvent.patient_mrn == mrn).delete()

    for event in events:
        db.add(models.ClinicalEvent(**event.model_dump(), patient_mrn=mrn))

    db.commit()

    # Return updated list
    return db.query(models.ClinicalEvent).filter(models.ClinicalEvent.patient_mrn == mrn).all()

@app.get("/api/patients/{mrn}/events", response_model=List[schemas.ClinicalEvent])
def read_events(mrn: str, db: Session = Depends(get_db), _=Depends(require_token)):
    events = db.query(models.ClinicalEvent).filter(models.ClinicalEvent.patient_mrn == mrn).all()
    return events

@app.post("/api/patients/{mrn}/outcomes", response_model=List[schemas.ClinicalOutcome])
def create_outcomes(mrn: str, outcomes: List[schemas.ClinicalOutcomeCreate], db: Session = Depends(get_db), _=Depends(require_token)):
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
def read_outcomes(mrn: str, db: Session = Depends(get_db), _=Depends(require_token)):
    outcomes = db.query(models.ClinicalOutcome).filter(models.ClinicalOutcome.patient_mrn == mrn).all()
    return outcomes

# Centre-wide calibration scalar (Center Calibration Report). Single row,
# id=1, so it applies to every device that loads the app rather than being
# stuck in whichever browser's localStorage clicked "Apply".
@app.get("/api/config/cl-scalar", response_model=schemas.CenterConfig)
def get_cl_scalar(db: Session = Depends(get_db), _=Depends(require_token)):
    cfg = db.query(models.CenterConfig).filter(models.CenterConfig.id == 1).first()
    if not cfg:
        return schemas.CenterConfig(cl_scalar=1.18, updated_at=None)
    return cfg

@app.post("/api/config/cl-scalar", response_model=schemas.CenterConfig)
def set_cl_scalar(payload: schemas.CenterConfigUpdate, db: Session = Depends(get_db), _=Depends(require_token)):
    cfg = db.query(models.CenterConfig).filter(models.CenterConfig.id == 1).first()
    if not cfg:
        cfg = models.CenterConfig(id=1, cl_scalar=payload.cl_scalar, updated_at=datetime.now(timezone.utc))
        db.add(cfg)
    else:
        cfg.cl_scalar = payload.cl_scalar
        cfg.updated_at = datetime.now(timezone.utc)
    db.commit()
    db.refresh(cfg)
    return cfg
