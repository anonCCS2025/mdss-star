#!/bin/bash
aws s3 cp /home/ubuntu/out.txt s3://sam-decoder-results-2025/$(date +%s).txt 
/home/ubuntu/cleanup.sh
