# Nif-Finder Web

Vercel-ready web interface for submitting protein FASTA sequences, running a
Nif-Finder-compatible analysis API, and visualizing predicted nif hits.

## Local development

```bash
npm install
npm run dev
```

Open http://localhost:3000.

To run the real Python/HMMER pipeline from the local web app:

```bash
NIF_FINDER_LOCAL_RUNNER=1 \
NIF_FINDER_ROOT=/path/to/Nif_finder/generl_bacteria \
NIF_FINDER_DB=/path/to/Nif_finder/generl_bacteria \
npm run dev
```

## Compute API

Set `NIF_FINDER_API_URL` in Vercel to connect the web app to a compute service
that runs the Python/HMMER Nif-Finder pipeline. The Vercel app posts:

```json
{
  "fasta": ">seq1\nM...",
  "jobs": 3,
  "cpu": 6
}
```

The compute service should return:

```json
{
  "records": [
    {
      "query": "seq1",
      "logEvalue": 164.62,
      "alignLength": 283,
      "queryLength": 296,
      "prediction": "nifH",
      "completeness": "Full"
    }
  ]
}
```

If `NIF_FINDER_API_URL` is not configured, the app returns a clear setup message
instead of pretending to run the HMMER pipeline.

Vercel is best used here for the web UI and visualization layer. Production
HMMER execution should run in a compute environment where the `hmmscan` binary,
database files, memory, and request duration can be controlled reliably.

This repository includes a Docker/FastAPI compute service in `../compute`.
Deploy that service first, then set `NIF_FINDER_API_URL` to its `/analyze`
endpoint and redeploy the Vercel app.
