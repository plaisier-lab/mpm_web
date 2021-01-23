#!/bin/bash
GLIOMA_SETTINGS=/home/me/Research/mpm_web/settings.cfg
PYTHONUNBUFFERED=0
docker-compose build app
docker-compose up --force-recreate --renew-anon-volumes
