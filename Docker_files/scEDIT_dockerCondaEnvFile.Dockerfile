# Use the official Python 3.9 image as the base
FROM ubuntu:22.04

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
    bzip2 \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
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

    # Download and install Anaconda
#RUN wget https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh -O /tmp/anaconda.sh && \
#ENV ANACONDA_VERSION=2023.03
ENV ANACONDA_DIR=/opt/anaconda3
COPY Anaconda3-2024.02-1-Linux-x86_64.sh /tmp/anaconda.sh 
RUN bash /tmp/anaconda.sh -b -p ${ANACONDA_DIR} && \
    rm /tmp/anaconda.sh && \
    ${ANACONDA_DIR}/bin/conda init && \
    ln -s ${ANACONDA_DIR}/bin/conda /usr/bin/conda
# Update PATH to include Conda
ENV PATH="${ANACONDA_DIR}/bin:${PATH}"

# Create a symbolic link for Python3 to Python3.9
#RUN ln -s /usr/bin/python3.9 /usr/bin/python3

# Upgrade pip for Python 3.9
#RUN python3.9 -m pip install --upgrade pip

# Install required Python packages
#COPY python_requirements.txt /docHome/python_requirements.txt
#RUN pip install --user --no-cache-dir -r /docHome/python_requirements.txt  && rm /docHome/python_requirements.txt 

# Set working directory
WORKDIR /docHome
# Copy the environment.yml file into the container
COPY scEDIT_condaEnv.yml /docHome/scEDIT_condaEnv.yml

# Update Conda and install the Conda environment
#RUN conda update -n base -c defaults conda -y 
RUN conda config --add channels defaults
RUN conda config --add channels bioconda 
RUN conda config --add channels conda-forge  
RUN conda env create -f /docHome/scEDIT_condaEnv.yml
RUN conda env update -n scEDIT -f /docHome/scEDIT_condaEnv.yml
# Activate the Conda environment by default
RUN echo "source ${ANACONDA_DIR}/bin/activate" >> ~/.bashrc
RUN echo "conda activate $(head -n 1 /docHome/scEDIT_condaEnv.yml | cut -d' ' -f2)" >> ~/.bashrc

#create user for scEDIT
#RUN useradd -m -s /bin/bash scEDIT && echo "scEDIT:scEDIT_1234" | chpasswd && adduser scEDIT sudo 
# Install required Python packages
#COPY python_requirements.txt /tmp/python_requirements.txt


# Set ownership of the working directory
#RUN chown -R scEDIT:scEDIT /docHome


# Clone GitHub repository
#RUN git clone https://github.com/GSbioinfo/scEDIT.git /docHome

# Set the default shell to Bash
SHELL ["/bin/bash", "--login", "-c"]

# Default command
CMD ["/bin/bash"]

# docker build --no-cache -t interactive_cond_scedit:1.0 --file Github_repos/scEDIT/Docker_files/scEDIT_dockerCondaEnvFile.Dockerfile /scratch/DockerImages/ --progress=plain
# docker save -o interactive_cond_scedit.tar interactive_cond_scedit:1.0
# docker load -i interactive_cond_scedit.tar