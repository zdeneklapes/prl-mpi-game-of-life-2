FROM ubuntu:20.04

# Install system packages
RUN DEBIAN_FRONTEND="noninteractive" apt-get update && apt-get -y install tzdata
RUN <<EOF
    apt-get update
    apt-get upgrade -y
    apt-get install -y \
        openssh-server  \
        sudo \
        ssh \
        python3-pip \
        python3-dev \
        python3-venv \
        libevent-dev \
        vim \
        gcc \
        g++ \
        gdb \
        clang \
        cmake \
        make \
        build-essential \
        autoconf \
        automake \
        valgrind \
        software-properties-common \
        libomp-dev \
        ninja-build \
        neovim \
        swig \
        locales-all \
        dos2unix \
        doxygen \
        fish \
        rsync \
        tar \
        tree \
        wget \
        mpich \
        bat
    apt-get clean
    rm -rf /var/lib/apt/lists/*
EOF

RUN <<EOF
    pip3 install --upgrade pip numpy matplotlib scipy pandas seaborn
EOF

WORKDIR /app
