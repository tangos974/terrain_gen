FROM ubuntu:latest

#https://github.com/EarthByte/PlateTectonicTools?tab=readme-ov-file

# Install packages for conda
RUN apt-get update && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git mercurial subversion libgl1 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and install Miniconda for python 3.10
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_24.3.0-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    /opt/conda/bin/conda clean -afy

# Set path to conda
ENV PATH /opt/conda/bin:$PATH

# Create and configure the pygplates environment and install additional tools
RUN /bin/bash -c "source /opt/conda/bin/activate && \
    conda create -n pygplates -c conda-forge python=3.10 pygplates cartopy platetectonictools -y" && \
    echo "source /opt/conda/bin/activate pygplates" >> ~/.bashrc

# Make RUN commands use the new environment:
SHELL ["/bin/bash", "--login", "-c"]