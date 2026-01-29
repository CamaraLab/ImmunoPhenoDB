#!/bin/sh
exec gunicorn --config /app/gunicorn_config.py app.wsgi:app