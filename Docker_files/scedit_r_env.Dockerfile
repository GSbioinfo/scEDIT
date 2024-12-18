# Use the Rocker R 4.4.1 image as the base
FROM rocker/r-ver:4.4.1

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC

# Update and install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libpng-dev \
    libjpeg-dev \
    zlib1g-dev \
    libudunits2-dev \
    libgdal-dev \
    libproj-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy R package installation script into the container
COPY Rpackages.R /tmp/Rpackages.R

# Install R packages
RUN Rscript /tmp/Rpackages.R && rm /tmp/Rpackages.R

# Set default working directory
WORKDIR /dockerHome

# Default command to start R
#CMD ["R"]
# Default interactive shell command
CMD ["/bin/bash"]

# docker build -t interactive_scedit_r_env:1.0 --file scEDIT/Docker_files/scedit_r_env.Dockerfile .
# docker save -o interactive_scedit_r_env.tar scedit_r_env:1.0
# docker load -i interactive_scedit_r_env.tar
# docker run -it -v /home/gajendra:/dockerHome interactive_scedit_r_env:1.0