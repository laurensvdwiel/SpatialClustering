FROM python:latest

WORKDIR /app

RUN pip install --upgrade pip setuptools wheel

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
COPY requirements.txt /app/
RUN pip install -r /app/requirements.txt