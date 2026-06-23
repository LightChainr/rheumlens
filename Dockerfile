FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /opt/rheumlens

RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential git \
    && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md LICENSE ./
COPY src ./src
COPY tests ./tests
COPY scripts ./scripts
COPY figures ./figures
COPY supplementary_tables ./supplementary_tables
COPY results_manifest ./results_manifest

RUN pip install --upgrade pip \
    && pip install -e ".[dev,io]"

CMD ["bash", "scripts/reproduce_minimal.sh"]
