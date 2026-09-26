@echo off
REM Run from the pipeline package root after editing pipeline_config.json
python run_full_pipeline.py --config pipeline_config.json --stages all
pause
