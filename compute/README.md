# Nif-Finder Compute API

Docker/FastAPI service for running the Python/HMMER Nif-Finder pipeline from the
Vercel web interface.

## Local run

```bash
cd /path/to/Nif_finder
python -m pip install -r compute/requirements.txt

NIF_FINDER_ROOT=$PWD/generl_bacteria \
NIF_FINDER_DB=$PWD/generl_bacteria \
uvicorn compute.app.main:app --host 0.0.0.0 --port 8000
```

Test:

```bash
curl http://localhost:8000/health
```

## Docker run

Build from the repository root so the Dockerfile can copy `generl_bacteria/`:

```bash
docker build -f compute/Dockerfile -t nif-finder-compute .
docker run --rm -p 8000:8000 nif-finder-compute
```

## Connect Vercel

After deploying this container to a compute host, set the Vercel web app
environment variable to the analyze endpoint:

```bash
cd web
npx vercel env add NIF_FINDER_API_URL production
# enter: https://your-compute-host.example.com/analyze
npx vercel --prod
```

For browser access, set `ALLOWED_ORIGINS` on the compute host to the Vercel app
origin, for example `https://web-theta-black-17.vercel.app`. During early
testing, the default `*` is permissive.
