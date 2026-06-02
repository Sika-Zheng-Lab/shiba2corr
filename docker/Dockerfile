FROM python:3.11-slim

LABEL maintainer="Naoto Kubota" \
      version="v1.0.0" \
      description="shiba2corr: weighted dPSI correlation analysis for Shiba"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/shiba2corr
COPY . /opt/shiba2corr

RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir .

WORKDIR /work

ENTRYPOINT ["shiba2corr"]
CMD ["--help"]
