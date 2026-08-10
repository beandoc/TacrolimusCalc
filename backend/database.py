from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import os

# Read from environment (Vercel) or fallback to local SQLite
SQLALCHEMY_DATABASE_URL = os.getenv("DATABASE_URL", "sqlite:///./tacrolimus.db")

# Fail fast rather than silently losing patient data.
#
# On a serverless platform the filesystem is ephemeral and per-invocation, so a
# SQLite file "works" — writes report success — and then vanishes, or is simply
# never visible to the next request. The frontend cannot tell this apart from a
# healthy backend: it gets 200s, shows "Saved to Database", and the record is
# gone. Refusing to start is the only honest failure mode.
_IS_SERVERLESS = bool(os.getenv("VERCEL") or os.getenv("AWS_LAMBDA_FUNCTION_NAME"))
if _IS_SERVERLESS and SQLALCHEMY_DATABASE_URL.startswith("sqlite"):
    raise RuntimeError(
        "DATABASE_URL is not set, so the API would fall back to SQLite on a "
        "serverless filesystem — writes would appear to succeed and then be lost. "
        "Set DATABASE_URL to a managed Postgres instance."
    )

# Only use check_same_thread for SQLite
connect_args = {}
if SQLALCHEMY_DATABASE_URL.startswith("sqlite"):
    connect_args = {"check_same_thread": False}

engine = create_engine(
    SQLALCHEMY_DATABASE_URL, connect_args=connect_args
)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()
