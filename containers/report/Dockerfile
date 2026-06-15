# Container for the LPS typing HTML report step (process: html_report).
# Pure Python; renders bin/generate_lps_report.py output. Build and push once,
# then pin the tag in nextflow.config (withName: html_report).
#
#   docker build -t ghcr.io/julianzaugg/lps_typing_report:0.1 containers/report
#   docker push ghcr.io/julianzaugg/lps_typing_report:0.1
#

FROM python:3.12-slim

LABEL org.opencontainers.image.source="https://github.com/julianzaugg/LPS_typing_Illumina"

RUN pip install --no-cache-dir \
        pandas==2.2.3 \
        matplotlib==3.9.2 \
        jinja2==3.1.4

# Matplotlib uses a non-interactive backend; ensure a writable config dir at runtime.
ENV MPLCONFIGDIR=/tmp/matplotlib

CMD ["python3"]
