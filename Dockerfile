FROM python:3.8-slim
WORKDIR /app

COPY requirements.txt .
RUN pip3 install -r requirements.txt

COPY s2p_setup.py .