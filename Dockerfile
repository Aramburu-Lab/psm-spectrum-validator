# The export script with the packages it was tested with, for a computer with Docker and no Python.
# The base is pinned by digest, so every build on one architecture starts from the same image.
FROM python:3.12-slim@sha256:b18b261e1b4974b4f321451ad5251ded85abb22c2b8b262da527dda1d031ab30
COPY requirements.txt constraints.txt /app/
RUN pip install --no-cache-dir --root-user-action=ignore -r /app/requirements.txt -c /app/constraints.txt
COPY scripts/python/ /app/scripts/python/
# The folder mounted here is where relative paths on the command line point.
WORKDIR /data
ENTRYPOINT ["python", "-I", "/app/scripts/python/export_mgf.py"]
