import base64
import hashlib
import hmac
import os
import subprocess
import tempfile
import threading
import time
from io import StringIO
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
from Bio import SeqIO
from fastapi import FastAPI, Header, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pygenomeviz import GenomeViz
from pydantic import BaseModel, Field


NIF_FINDER_ROOT = Path(os.environ.get("NIF_FINDER_ROOT", "/app/generl_bacteria"))
NIF_FINDER_DB = Path(os.environ.get("NIF_FINDER_DB", str(NIF_FINDER_ROOT)))
NIF_FINDER_SCRIPT = NIF_FINDER_ROOT / "Nif_finderv0_24.py"
PYTHON_BIN = os.environ.get("NIF_FINDER_PYTHON", "python")
MAX_FASTA_BYTES = int(os.environ.get("MAX_FASTA_BYTES", str(10 * 1024 * 1024)))
MAX_GENBANK_BYTES = int(os.environ.get("MAX_GENBANK_BYTES", str(30 * 1024 * 1024)))
REQUEST_TIMEOUT_SECONDS = max(600, int(os.environ.get("REQUEST_TIMEOUT_SECONDS", "600")))
MAX_CONCURRENT_ANALYSES = int(os.environ.get("MAX_CONCURRENT_ANALYSES", "2"))
QUEUE_WAIT_SECONDS = int(os.environ.get("QUEUE_WAIT_SECONDS", "240"))
NIF_FINDER_API_KEY = os.environ.get("NIF_FINDER_API_KEY")
TOKEN_MAX_AGE_SECONDS = int(os.environ.get("TOKEN_MAX_AGE_SECONDS", "600"))
NIF_GENES = {"nifH", "nifD", "nifK", "nifE", "nifN", "nifB"}
GENE_COLORS = {
    "nifH": "#E74C3C",
    "nifD": "#3498DB",
    "nifK": "#2ECC71",
    "nifE": "#F39C12",
    "nifN": "#9B59B6",
    "nifB": "#1ABC9C",
}
LOCAL_CONTEXT_PADDING = int(os.environ.get("LOCAL_CONTEXT_PADDING", "10000"))
LOCAL_CONTEXT_MERGE_DISTANCE = int(os.environ.get("LOCAL_CONTEXT_MERGE_DISTANCE", "20000"))
OVERVIEW_HIDE_COVERAGE = float(os.environ.get("OVERVIEW_HIDE_COVERAGE", "0.5"))
TRACK_LINE_KWS = {"color": "#20242a", "lw": 0.45, "zorder": 0}


class AnalyzeRequest(BaseModel):
    fasta: str = Field(..., min_length=1)
    genbank: str | None = None
    jobs: int = Field(default=1, ge=1, le=4)
    cpu: int = Field(default=4, ge=1, le=12)
    plot: bool = True
    evalue: float = Field(default=1e-10, gt=0)


class ResultRecord(BaseModel):
    query: str
    logEvalue: float
    alignLength: int
    queryLength: int | None
    prediction: str
    completeness: str
    operonLabel: str | None = None


class AnalyzeResponse(BaseModel):
    records: list[ResultRecord]
    plotPngBase64: str | None = None
    plotMessage: str | None = None
    genomicContextOverviewSvg: str | None = None
    genomicContextLocalSvg: str | None = None
    genomicContextMessage: str | None = None
    runner: str = "nif-finder-compute"


class CdsFeature(BaseModel):
    contig: str
    contig_length: int
    cds_index: int
    start: int
    end: int
    strand: int
    protein_id: str | None = None
    locus_tag: str | None = None
    old_locus_tag: str | None = None
    gene: str | None = None
    translation: str | None = None


class MatchedHit(BaseModel):
    feature: CdsFeature
    prediction: str
    query: str
    completeness: str
    operon_label: str | None = None


app = FastAPI(title="Nif-Finder Compute API", version="0.1.0")
analysis_semaphore = threading.BoundedSemaphore(MAX_CONCURRENT_ANALYSES)

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
        expected_with_operon = [*expected, "Operon"]
        if header not in (expected, expected_with_operon):
            raise ValueError(f"Unexpected Nif-Finder output header: {header}")

        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            query, log_evalue, align_len, query_len, prediction, completeness = fields[:6]
            operon_label = fields[6] if len(fields) > 6 else ""
            records.append(
                {
                    "query": query,
                    "logEvalue": float(log_evalue),
                    "alignLength": int(align_len),
                    "queryLength": None if query_len == "N/A" else int(query_len),
                    "prediction": prediction,
                    "completeness": completeness,
                    "operonLabel": operon_label or None,
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


def validate_genbank(genbank: str | None) -> str | None:
    if genbank is None or not genbank.strip():
        return None
    genbank = genbank.strip()
    if len(genbank.encode("utf-8")) > MAX_GENBANK_BYTES:
        raise HTTPException(status_code=413, detail=f"GenBank input exceeds {MAX_GENBANK_BYTES} bytes.")
    return genbank + "\n"


def normalize_identifier(value: str | None) -> str:
    parts = (value or "").strip().lstrip(">").split()
    return parts[0] if parts else ""


def query_identifier_candidates(query: str) -> list[str]:
    first = normalize_identifier(query)
    pieces = [first]
    for separator in ("|", ";", ","):
        if separator in first:
            pieces.extend(part.strip() for part in first.split(separator))
    if "_" in first:
        pieces.append(first.rsplit("_", 1)[0])
    seen: set[str] = set()
    candidates: list[str] = []
    for piece in pieces:
        if piece and piece not in seen and piece not in NIF_GENES:
            candidates.append(piece)
            seen.add(piece)
    return candidates


def query_cds_index_candidate(query: str) -> int | None:
    first = normalize_identifier(query)
    suffix = first.rsplit("_", 1)[-1]
    if not suffix.isdigit():
        return None
    index = int(suffix)
    return index if index > 0 else None


def qualifier_value(feature: Any, key: str) -> str | None:
    values = feature.qualifiers.get(key, [])
    if not values:
        return None
    return str(values[0])


def parse_fasta_sequences(fasta: str) -> dict[str, str]:
    sequences: dict[str, str] = {}
    for record in SeqIO.parse(StringIO(fasta), "fasta"):
        sequences[normalize_identifier(record.id)] = str(record.seq).replace("*", "").upper()
    return sequences


def parse_genbank_cds_features(genbank: str) -> list[CdsFeature]:
    features: list[CdsFeature] = []
    try:
        records = list(SeqIO.parse(StringIO(genbank), "genbank"))
    except Exception as exc:
        raise ValueError(f"Could not parse GenBank file: {exc}") from exc
    if not records:
        raise ValueError("No GenBank records were found.")

    for record in records:
        contig = record.id or record.name or "record"
        contig_length = len(record.seq)
        if contig_length <= 0:
            continue
        cds_index = 0
        for feature in record.features:
            if feature.type != "CDS" or feature.location is None:
                continue
            cds_index += 1
            start = int(feature.location.start)
            end = int(feature.location.end)
            if end <= start:
                continue
            strand = feature.location.strand if feature.location.strand in (-1, 1) else 1
            features.append(
                CdsFeature(
                    contig=contig,
                    contig_length=contig_length,
                    cds_index=cds_index,
                    start=max(0, start),
                    end=min(contig_length, end),
                    strand=strand,
                    protein_id=qualifier_value(feature, "protein_id"),
                    locus_tag=qualifier_value(feature, "locus_tag"),
                    old_locus_tag=qualifier_value(feature, "old_locus_tag"),
                    gene=qualifier_value(feature, "gene"),
                    translation=(qualifier_value(feature, "translation") or "").replace("*", "").upper() or None,
                )
            )
    return features


def match_nif_hits_to_genbank(
    records: list[dict[str, Any]],
    features: list[CdsFeature],
    fasta_sequences: dict[str, str] | None = None,
) -> list[MatchedHit]:
    by_protein_id: dict[str, list[CdsFeature]] = {}
    by_locus_tag: dict[str, list[CdsFeature]] = {}
    by_old_locus_tag: dict[str, list[CdsFeature]] = {}
    by_gene: dict[str, list[CdsFeature]] = {}
    by_cds_index: dict[int, CdsFeature] = {}
    by_translation: dict[str, list[CdsFeature]] = {}
    for feature in features:
        by_cds_index.setdefault(feature.cds_index, feature)
        if feature.translation:
            by_translation.setdefault(feature.translation, []).append(feature)
        for mapping, value in (
            (by_protein_id, feature.protein_id),
            (by_locus_tag, feature.locus_tag),
            (by_old_locus_tag, feature.old_locus_tag),
            (by_gene, feature.gene),
        ):
            key = normalize_identifier(value)
            if key:
                mapping.setdefault(key, []).append(feature)

    matches: list[MatchedHit] = []
    used: set[tuple[str, int, int, str]] = set()
    prediction_counts: dict[str, int] = {}
    for record in records:
        prediction = str(record.get("prediction") or "")
        if prediction in NIF_GENES:
            prediction_counts[prediction] = prediction_counts.get(prediction, 0) + 1

    for record in records:
        prediction = str(record.get("prediction") or "")
        if prediction not in NIF_GENES:
            continue
        query = str(record.get("query") or "")
        candidates = query_identifier_candidates(query)
        found: list[CdsFeature] = []
        for mapping in (by_protein_id, by_locus_tag, by_old_locus_tag):
            for candidate in candidates:
                found = mapping.get(candidate, [])
                if found:
                    break
            if found:
                break
        if not found:
            exact_query = normalize_identifier(query)
            gene_matches = by_gene.get(exact_query, [])
            if len(gene_matches) == 1:
                found = gene_matches
            elif exact_query == prediction and len(by_gene.get(prediction, [])) == 1:
                found = by_gene[prediction]
            elif prediction_counts.get(prediction) == 1 and len(by_gene.get(prediction, [])) == 1:
                found = by_gene[prediction]
        if not found:
            cds_index = query_cds_index_candidate(query)
            feature = by_cds_index.get(cds_index) if cds_index else None
            if feature and normalize_identifier(feature.gene) == prediction:
                found = [feature]
        if not found and fasta_sequences:
            for candidate in candidates:
                sequence = fasta_sequences.get(candidate)
                if not sequence:
                    continue
                translated_matches = [
                    feature
                    for feature in by_translation.get(sequence, [])
                    if normalize_identifier(feature.gene) in ("", prediction)
                ]
                if translated_matches:
                    found = translated_matches
                    break

        for feature in found:
            key = (feature.contig, feature.start, feature.end, prediction)
            if key in used:
                continue
            used.add(key)
            matches.append(
                MatchedHit(
                    feature=feature,
                    prediction=prediction,
                    query=query,
                    completeness=str(record.get("completeness") or ""),
                    operon_label=record.get("operonLabel") or None,
                )
            )
    return matches


def svg_to_text(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    svg_start = text.find("<svg")
    return text[svg_start:] if svg_start >= 0 else text


def add_feature(track: Any, feature: CdsFeature, *, label: str = "", color: str = "#c7cdd4") -> None:
    track.add_feature(
        feature.start,
        feature.end,
        feature.strand,
        plotstyle="bigarrow",
        label=label,
        fc=color,
        ec=color,
        lw=0,
        zorder=2,
        text_kws={"size": 9, "rotation": 0, "vpos": "top", "hpos": "center"} if label else None,
    )


def local_hit_suffix(hit: MatchedHit) -> str:
    if hit.completeness == "Full_operon":
        return "(operon)"
    if hit.completeness == "Full":
        return "(full)"
    if hit.completeness == "Fragment":
        return "(frag.)"
    return ""


def hit_label(hit: MatchedHit) -> str:
    if hit.completeness == "Full_operon" and hit.operon_label:
        return hit.operon_label
    return hit.prediction


def overview_hit_label(hit: MatchedHit) -> str:
    label = hit_label(hit)
    if label.startswith("nif"):
        return label.replace("nif", "")
    return label


def add_overview_position_labels(track: Any, contig_length: int) -> None:
    track.add_text(0, "1", size=8, vpos="bottom", hpos="left", ymargin=0.25, rotation=0, color="#20242a")
    track.add_text(
        contig_length,
        f"{contig_length:,}",
        size=8,
        vpos="bottom",
        hpos="right",
        ymargin=0.25,
        rotation=0,
        color="#20242a",
    )


def add_hit_suffix(track: Any, feature: CdsFeature, hit: MatchedHit) -> None:
    suffix = local_hit_suffix(hit)
    if not suffix:
        return
    track.add_text(
        (feature.start + feature.end) / 2,
        suffix,
        size=8,
        vpos="bottom",
        hpos="center",
        ymargin=0.2,
        rotation=0,
        color=GENE_COLORS.get(hit.prediction, "#0f766e"),
    )


def clip_feature_to_region(feature: CdsFeature, start: int, end: int) -> CdsFeature:
    return feature.model_copy(update={"start": max(feature.start, start), "end": min(feature.end, end)})


def build_overview_svg(matches: list[MatchedHit], out_path: Path) -> None:
    contigs = sorted(
        {
            match.feature.contig: match.feature.contig_length
            for match in matches
        }.items(),
        key=lambda item: item[0],
    )
    fig_height = max(0.45, min(0.8, 3.8 / max(1, len(contigs))))
    gv = GenomeViz(fig_width=12, fig_track_height=fig_height, show_axis=False, theme="light")
    for contig, contig_length in contigs:
        track = gv.add_feature_track(contig, contig_length, labelsize=10, line_kws=TRACK_LINE_KWS)
        add_overview_position_labels(track, contig_length)
        for match in sorted((m for m in matches if m.feature.contig == contig), key=lambda m: m.feature.start):
            add_feature(
                track,
                match.feature,
                label=overview_hit_label(match),
                color=GENE_COLORS.get(match.prediction, "#0f766e"),
            )
    gv.savefig(out_path, dpi=140, pad_inches=0.25)


def group_local_regions(matches: list[MatchedHit]) -> list[tuple[str, int, int, list[MatchedHit]]]:
    regions: list[tuple[str, int, int, list[MatchedHit]]] = []
    by_contig: dict[str, list[MatchedHit]] = {}
    for match in matches:
        by_contig.setdefault(match.feature.contig, []).append(match)

    for contig, contig_matches in sorted(by_contig.items()):
        ordered = sorted(contig_matches, key=lambda match: match.feature.start)
        current: list[MatchedHit] = []
        current_end = -1
        for match in ordered:
            if not current or match.feature.start - current_end <= LOCAL_CONTEXT_MERGE_DISTANCE:
                current.append(match)
                current_end = max(current_end, match.feature.end)
            else:
                contig_length = current[0].feature.contig_length
                start = max(0, min(m.feature.start for m in current) - LOCAL_CONTEXT_PADDING)
                end = min(contig_length, max(m.feature.end for m in current) + LOCAL_CONTEXT_PADDING)
                regions.append((contig, start, end, current))
                current = [match]
                current_end = match.feature.end
        if current:
            contig_length = current[0].feature.contig_length
            start = max(0, min(m.feature.start for m in current) - LOCAL_CONTEXT_PADDING)
            end = min(contig_length, max(m.feature.end for m in current) + LOCAL_CONTEXT_PADDING)
            regions.append((contig, start, end, current))
    return regions


def should_hide_whole_genome_view(matches: list[MatchedHit]) -> bool:
    regions = group_local_regions(matches)
    by_contig: dict[str, list[tuple[str, int, int, list[MatchedHit]]]] = {}
    for region in regions:
        by_contig.setdefault(region[0], []).append(region)

    if not by_contig:
        return False

    for contig_regions in by_contig.values():
        if len(contig_regions) != 1:
            return False
        _contig, start, end, region_matches = contig_regions[0]
        contig_length = region_matches[0].feature.contig_length
        if contig_length <= 0:
            return False
        coverage = (end - start) / contig_length
        if coverage < OVERVIEW_HIDE_COVERAGE:
            return False
    return True


def build_local_context_svg(matches: list[MatchedHit], features: list[CdsFeature], out_path: Path) -> None:
    regions = group_local_regions(matches)
    fig_height = max(0.55, min(0.9, 4.5 / max(1, len(regions))))
    gv = GenomeViz(fig_width=12, fig_track_height=fig_height, show_axis=False, theme="light")
    hit_keys = {
        (match.feature.contig, match.feature.start, match.feature.end): match
        for match in matches
    }
    contig_region_totals = {contig: sum(1 for region in regions if region[0] == contig) for contig, *_ in regions}
    contig_region_counts: dict[str, int] = {}
    for contig, start, end, _region_matches in regions:
        contig_region_counts[contig] = contig_region_counts.get(contig, 0) + 1
        region_name = f"{contig}: {start + 1:,}-{end:,}"
        if contig_region_totals[contig] > 1:
            region_name = f"{contig} region {contig_region_counts[contig]}: {start + 1:,}-{end:,}"
        track = gv.add_feature_track(region_name, (start, end), labelsize=10, line_kws=TRACK_LINE_KWS)
        overlapping = [
            feature
            for feature in features
            if feature.contig == contig and feature.end >= start and feature.start <= end
        ]
        for feature in sorted(overlapping, key=lambda item: item.start):
            hit = hit_keys.get((feature.contig, feature.start, feature.end))
            clipped_feature = clip_feature_to_region(feature, start, end)
            if hit:
                add_feature(
                    track,
                    clipped_feature,
                    label=hit_label(hit),
                    color=GENE_COLORS.get(hit.prediction, "#0f766e"),
                )
                add_hit_suffix(track, clipped_feature, hit)
            else:
                add_feature(track, clipped_feature)
    gv.savefig(out_path, dpi=140, pad_inches=0.25)


def build_genomic_context(
    genbank: str | None,
    records: list[dict[str, Any]],
    tmp_dir: Path,
    fasta: str | None = None,
) -> dict[str, str | None]:
    if not genbank:
        return {}
    try:
        features = parse_genbank_cds_features(genbank)
        fasta_sequences = parse_fasta_sequences(fasta) if fasta else None
        matches = match_nif_hits_to_genbank(records, features, fasta_sequences)
        if not matches:
            return {
                "genomicContextMessage": (
                    "A GenBank file was provided, but no Nif-Finder hits could be matched to CDS "
                    "protein_id, locus_tag, old_locus_tag, a unique gene qualifier, or CDS translation."
                )
            }
        overview_path = tmp_dir / "genomic_context_overview.svg"
        local_path = tmp_dir / "genomic_context_local.svg"
        hide_overview = should_hide_whole_genome_view(matches)
        if not hide_overview:
            build_overview_svg(matches, overview_path)
        build_local_context_svg(matches, features, local_path)
        return {
            "genomicContextOverviewSvg": None if hide_overview else svg_to_text(overview_path),
            "genomicContextLocalSvg": svg_to_text(local_path),
            "genomicContextMessage": None,
        }
    except Exception as exc:
        return {"genomicContextMessage": f"Genomic context visualization failed: {exc}"}


def validate_analysis_token(token: str | None) -> bool:
    if not token or not NIF_FINDER_API_KEY:
        return False
    parts = token.split(".")
    if len(parts) != 3:
        return False
    expires, nonce, signature = parts
    if not expires.isdigit():
        return False
    expires_at = int(expires)
    now = int(time.time())
    if expires_at < now or expires_at > now + TOKEN_MAX_AGE_SECONDS:
        return False
    message = f"{expires}.{nonce}".encode("utf-8")
    expected = hmac.new(NIF_FINDER_API_KEY.encode("utf-8"), message, hashlib.sha256).hexdigest()
    return hmac.compare_digest(signature, expected)


@app.get("/health")
def health() -> dict[str, str | bool | int]:
    hmmscan = subprocess.run(["hmmscan", "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    status = "ok" if hmmscan.returncode == 0 and NIF_FINDER_SCRIPT.is_file() else "not-ready"
    try:
        import matplotlib  # noqa: F401

        matplotlib_available = True
    except ImportError:
        matplotlib_available = False

    return {
        "status": status,
        "nifFinderRoot": str(NIF_FINDER_ROOT),
        "nifFinderDb": str(NIF_FINDER_DB),
        "matplotlib": matplotlib_available,
        "requestTimeoutSeconds": REQUEST_TIMEOUT_SECONDS,
    }


@app.post("/analyze", response_model=AnalyzeResponse)
def analyze(
    request: AnalyzeRequest,
    x_api_key: str | None = Header(default=None),
    x_analysis_token: str | None = Header(default=None),
) -> AnalyzeResponse:
    api_key_valid = bool(NIF_FINDER_API_KEY and hmac.compare_digest(x_api_key or "", NIF_FINDER_API_KEY))
    token_valid = validate_analysis_token(x_analysis_token)
    if NIF_FINDER_API_KEY and not (api_key_valid or token_valid):
        raise HTTPException(status_code=401, detail="Invalid or missing API key.")

    fasta = validate_fasta(request.fasta)
    genbank = validate_genbank(request.genbank)

    if not NIF_FINDER_SCRIPT.is_file():
        raise HTTPException(status_code=500, detail=f"Nif-Finder script not found: {NIF_FINDER_SCRIPT}")
    if not NIF_FINDER_DB.is_dir():
        raise HTTPException(status_code=500, detail=f"Nif-Finder DB directory not found: {NIF_FINDER_DB}")

    acquired = analysis_semaphore.acquire(timeout=QUEUE_WAIT_SECONDS)
    if not acquired:
        raise HTTPException(
            status_code=429,
            detail="The compute server is busy. Please try again in a few minutes.",
        )

    try:
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
                "--evalue",
                str(request.evalue),
            ]
            if request.plot:
                command.append("-p")
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

            plot_path = Path(f"{out_prefix}_scatter.png")
            plot_png_base64 = None
            plot_message = None
            if plot_path.is_file():
                plot_png_base64 = base64.b64encode(plot_path.read_bytes()).decode("ascii")
            elif not request.plot:
                plot_message = "Plot output was disabled for this run."
            else:
                output = "\n".join(part for part in [completed.stdout.strip(), completed.stderr.strip()] if part)
                plot_message = output[-1000:] if output else "Nif-Finder did not produce a scatter plot."

            genomic_context = build_genomic_context(genbank, records, tmp_dir, fasta)
    finally:
        analysis_semaphore.release()

    return AnalyzeResponse(
        records=records,
        plotPngBase64=plot_png_base64,
        plotMessage=plot_message,
        **genomic_context,
    )
