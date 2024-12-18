# Use the official Python 3.9 image as the base
FROM python:3.9-slim

# Set environment variables to avoid prompts during build
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1

# Install system dependencies
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
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    build-essential \
    libssl-dev \
    libffi-dev \
    wget \
    curl \
    python3-venv \
    python3-distutils \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Create a symbolic link for Python3 to Python3.9
#RUN ln -s /usr/bin/python3.9 /usr/bin/python3

# Upgrade pip for Python 3.9
RUN python3.9 -m pip install --upgrade pip

# Install required Python packages
COPY python_requirements.txt /docHome/python_requirements.txt
RUN pip install --user --no-cache-dir -r /docHome/python_requirements.txt  && rm /docHome/python_requirements.txt 

#create user for scEDIT
#RUN useradd -m -s /bin/bash scEDIT && echo "scEDIT:scEDIT_1234" | chpasswd && adduser scEDIT sudo 
# Install required Python packages
#COPY python_requirements.txt /tmp/python_requirements.txt

# Set working directory
WORKDIR /docHome
# Set ownership of the working directory
#RUN chown -R scEDIT:scEDIT /docHome



# Clone GitHub repository
#RUN git clone https://github.com/GSbioinfo/scEDIT.git /docHome

# Default command
CMD ["/bin/bash"]

# docker build --no-cache -t interactive_run_scedit:1.0 --file /scEDIT/Docker_files/scEDIT_dockerEnvFile.Dockerfile . --progress=plain 
# docker run -it -v /scratch/:/docHome/scratch/ interactive_run_scedit:1.0
