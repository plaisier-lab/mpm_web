FROM python:3.6

EXPOSE 5000

WORKDIR /app

COPY requirements.txt /app
RUN pip install -r requirements.txt

COPY app.py /app
COPY models.py /app
COPY populateDatabase.py /app
COPY templates /app/templates
COPY static /app/static

COPY data /app/data

CMD python app.py
