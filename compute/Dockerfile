FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    NIF_FINDER_ROOT=/app/generl_bacteria \
    NIF_FINDER_DB=/app/generl_bacteria

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends hmmer \
    && rm -rf /var/lib/apt/lists/*

COPY compute/requirements.txt /app/compute/requirements.txt
RUN pip install --no-cache-dir -r /app/compute/requirements.txt

COPY generl_bacteria /app/generl_bacteria
COPY compute/app /app/compute/app

EXPOSE 8000

CMD ["uvicorn", "compute.app.main:app", "--host", "0.0.0.0", "--port", "8000"]
