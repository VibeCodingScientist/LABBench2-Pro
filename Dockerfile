FROM python:3.11-slim

WORKDIR /app

# System deps for biopython/scipy
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc g++ libpq-dev && \
    rm -rf /var/lib/apt/lists/*

COPY pyproject.toml .
RUN pip install --no-cache-dir .

COPY . .

# Default: run the full pipeline
ENTRYPOINT ["python", "-m"]
CMD ["src.api"]
