# Use the latest Ubuntu LTS as the base image
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC

# Update and install basic dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libffi-dev \
    libgomp1 \
    libstdc++6 \
    libgcc-11-dev \
    libgfortran5 \
    libblas-dev \
    liblapack-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libopenblas-dev \
    libtbb-dev \
    ca-certificates \
    pigz \
    ncurses-dev \
    xz-utils \
    tk-dev \
    sqlite3 \
    libreadline-dev \
    libssl-dev \
    libgd-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Python 3.9 and pip
RUN add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y python3.9 python3.9-distutils python3.9-dev && \
    curl -sS https://bootstrap.pypa.io/get-pip.py | python3.9 && \
    ln -s /usr/bin/python3.9 /usr/bin/python3 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*



# Install required Python packages
COPY python_requirements.txt /tmp/python_requirements.txt
RUN pip install --no-cache-dir -r /tmp/python_requirements.txt && rm /tmp/python_requirements.txt


# Clone GitHub repository
RUN git clone https://github.com/your-repo/example-repo.git /opt/example-repo

# Set working directory
WORKDIR /opt/example-repo

# Default command
CMD ["/bin/bash"]