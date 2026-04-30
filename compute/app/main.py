import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field


NIF_FINDER_ROOT = Path(os.environ.get("NIF_FINDER_ROOT", "/app/generl_bacteria"))
NIF_FINDER_DB = Path(os.environ.get("NIF_FINDER_DB", str(NIF_FINDER_ROOT)))
NIF_FINDER_SCRIPT = NIF_FINDER_ROOT / "Nif_finderv0_24.py"
PYTHON_BIN = os.environ.get("NIF_FINDER_PYTHON", "python")
MAX_FASTA_BYTES = int(os.environ.get("MAX_FASTA_BYTES", str(10 * 1024 * 1024)))
REQUEST_TIMEOUT_SECONDS = int(os.environ.get("REQUEST_TIMEOUT_SECONDS", "120"))


class AnalyzeRequest(BaseModel):
    fasta: str = Field(..., min_length=1)
    jobs: int = Field(default=3, ge=1, le=6)
    cpu: int = Field(default=6, ge=1, le=32)


class ResultRecord(BaseModel):
    query: str
    logEvalue: float
    alignLength: int
    queryLength: int | None
    prediction: str
    completeness: str


class AnalyzeResponse(BaseModel):
    records: list[ResultRecord]
    runner: str = "nif-finder-compute"


app = FastAPI(title="Nif-Finder Compute API", version="0.1.0")

allowed_origins = [
    origin.strip()
    for origin in os.environ.get("ALLOWED_ORIGINS", "*").split(",")
    if origin.strip()
]
app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_credentials=False,
    allow_methods=["POST", "GET"],
    allow_headers=["*"],
)


def parse_results_tsv(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        expected = ["Query", "-log_Evalue", "Align_Len", "Query_Length", "Prediction", "Completeness"]
        if header != expected:
            raise ValueError(f"Unexpected Nif-Finder output header: {header}")

        for line in handle:
            if not line.strip():
                continue
            query, log_evalue, align_len, query_len, prediction, completeness = line.rstrip("\n").split("\t")
            records.append(
                {
                    "query": query,
                    "logEvalue": float(log_evalue),
                    "alignLength": int(align_len),
                    "queryLength": None if query_len == "N/A" else int(query_len),
                    "prediction": prediction,
                    "completeness": completeness,
                }
            )
    return records


def validate_fasta(fasta: str) -> str:
    fasta = fasta.strip()
    if not fasta:
        raise HTTPException(status_code=400, detail="Please provide a protein FASTA sequence.")
    if not fasta.startswith(">"):
        raise HTTPException(status_code=400, detail="Input must be FASTA formatted and start with '>'.")
    if len(fasta.encode("utf-8")) > MAX_FASTA_BYTES:
        raise HTTPException(status_code=413, detail=f"FASTA input exceeds {MAX_FASTA_BYTES} bytes.")
    return fasta + "\n"


@app.get("/health")
def health() -> dict[str, str]:
    hmmscan = subprocess.run(["hmmscan", "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    status = "ok" if hmmscan.returncode == 0 and NIF_FINDER_SCRIPT.is_file() else "not-ready"
    return {
        "status": status,
        "nifFinderRoot": str(NIF_FINDER_ROOT),
        "nifFinderDb": str(NIF_FINDER_DB),
    }


@app.post("/analyze", response_model=AnalyzeResponse)
def analyze(request: AnalyzeRequest) -> AnalyzeResponse:
    fasta = validate_fasta(request.fasta)

    if not NIF_FINDER_SCRIPT.is_file():
        raise HTTPException(status_code=500, detail=f"Nif-Finder script not found: {NIF_FINDER_SCRIPT}")
    if not NIF_FINDER_DB.is_dir():
        raise HTTPException(status_code=500, detail=f"Nif-Finder DB directory not found: {NIF_FINDER_DB}")

    with tempfile.TemporaryDirectory(prefix="nif-finder-api-") as tmp:
        tmp_dir = Path(tmp)
        query_path = tmp_dir / "query.faa"
        out_prefix = tmp_dir / "result"
        query_path.write_text(fasta)

        command = [
            PYTHON_BIN,
            str(NIF_FINDER_SCRIPT),
            "-q",
            str(query_path),
            "-o",
            str(out_prefix),
            "--jobs",
            str(request.jobs),
            "--cpu",
            str(request.cpu),
        ]
        env = {
            **os.environ,
            "NIF_FINDER_DB": str(NIF_FINDER_DB),
        }

        try:
            completed = subprocess.run(
                command,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=REQUEST_TIMEOUT_SECONDS,
                check=True,
            )
        except subprocess.TimeoutExpired as exc:
            raise HTTPException(status_code=504, detail=f"Nif-Finder timed out after {REQUEST_TIMEOUT_SECONDS}s.") from exc
        except subprocess.CalledProcessError as exc:
            detail = exc.stderr.strip() or exc.stdout.strip() or str(exc)
            raise HTTPException(status_code=500, detail=f"Nif-Finder failed: {detail}") from exc

        result_path = Path(f"{out_prefix}.txt")
        if not result_path.is_file():
            raise HTTPException(status_code=500, detail="Nif-Finder did not produce a result table.")

        try:
            records = parse_results_tsv(result_path)
        except Exception as exc:
            raise HTTPException(status_code=500, detail=f"Failed to parse Nif-Finder results: {exc}") from exc

    return AnalyzeResponse(records=records)
