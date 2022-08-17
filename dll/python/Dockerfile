FROM python:3.9-slim AS builder

ENV PYTHONUNBUFFERED 1

RUN apt-get update -y && \
    apt-get install -y automake build-essential libtool wget make cmake

RUN wget https://mirror.ibcp.fr/pub/gnu/gsl/gsl-2.7.tar.gz && \
    wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz && \
    wget https://github.com/stevengj/nlopt/archive/v2.7.1.tar.gz && \
    \
    tar -xvf \gsl-2.7.tar.gz && \
    cd gsl-2.7 && ./configure && make && make check && make install && cd / && \
    \
    tar -xvf eigen-3.4.0.tar.gz && \
    cd eigen-3.4.0 && mkdir build && cd build && cmake .. && make install && cd / && \
    \
    tar -xvf \v2.7.1.tar.gz && \
    cd nlopt-2.7.1/ && mkdir build && cd build && cmake .. && make install && \
    \
    rm -rf /eigen-3.4.0 /eigen-3.4.0.tar.gz /nlopt-2.7.1/ /v2.7.1.tar.gz

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

COPY . /app/

RUN cd /app && \
    mkdir -p ../BMD_DLL_COMPILE && \
    ./create_dll_compile.bash && \
    cd /BMD_DLL_COMPILE && \
    make -j4 && \
    make install && \
    python /app/dll/python/dll_test.py

# -------------------------------------------------------------------------------------------------
# copy outputs into clean docker build
FROM python:3.9-slim

COPY --from=builder /usr/local/lib/libDRBMD.* /usr/local/lib/
COPY --from=builder /usr/local/lib/libnlopt* /usr/local/lib/
COPY --from=builder /app/dll/python/dll_test.py /app/dll/python/dll_test.py

RUN groupadd -g 555 -r app && \
    useradd -u 555 -r -g app app && \
    chown -R app:app /app

# demo using non-root user

USER app

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

RUN ldd /usr/local/lib/libDRBMD.so && \
    python /app/dll/python/dll_test.py
