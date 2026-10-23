FROM rocker/r-ver:4.5.0
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-venv python3-pip build-essential gfortran libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev libtbb-dev libfontconfig1-dev libfreetype6-dev libpng-dev libjpeg-dev libtiff-dev zlib1g-dev && rm -rf /var/lib/apt/lists/*
WORKDIR /analysis
COPY . /analysis
RUN bash scripts/setup_linux.sh
CMD ["bash", "scripts/run_full.sh"]
