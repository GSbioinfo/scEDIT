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


# Install Perl and required Perl modules
RUN apt-get update && \
    apt-get install -y perl && \
    curl -L https://cpanmin.us | perl - App::cpanminus && \
    cpanm Math::Bezier Math::VecStat Data::Dumper List::MoreUtils GD::SVG Color::Rgb && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Clone GitHub repository
RUN git clone https://github.com/GSbioinfo/PetalPlot.git /opt/

# Set working directory
WORKDIR /opt/

# Default command
CMD ["/bin/bash"]