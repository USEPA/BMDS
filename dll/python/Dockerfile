FROM python:3.11-slim AS builder

ENV PYTHONUNBUFFERED 1

RUN apt-get update -y && \
    apt-get install -y automake build-essential libtool make cmake libgslcblas0 libgsl-dev

COPY ./vendor /vendor
RUN cd /vendor && \
    \
    tar -xvf nlopt-2.7.1.tar.gz && \
    tar -xvf eigen-3.4.0.tar.gz && \
    \
    cd /vendor/nlopt-2.7.1 && mkdir build && cd build && cmake .. && make install && \
    cd /vendor/eigen-3.4.0 && mkdir build && cd build && cmake .. && make install && \
    \
    rm -rf /vendor

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

COPY ./dll /app/dll

COPY ./create_dll_compile.bash.template /app/create_dll_compile.bash
COPY . /app
RUN cd /app && \
    mkdir -p ../BMD_DLL_COMPILE && \
    ./create_dll_compile.bash && \
    cd /BMD_DLL_COMPILE && \
    make -j$(nproc) && \
    strip --strip-all ./.libs/libDRBMD.so.0.0.0 && \
    make install && \
    python /app/dll/python/dll_test.py

# -----------------------------------------------------------------------------
# copy outputs into clean docker build
FROM python:3.11-slim

COPY --from=builder /usr/local/lib/libDRBMD.* /usr/local/lib/
COPY --from=builder /usr/local/lib/libnlopt* /usr/local/lib/
COPY --from=builder /app/dll/python/dll_test.py /app/dll/python/dll_test.py

RUN apt-get update -y && \
    apt-get install -y libgslcblas0 libgsl-dev

RUN groupadd -g 555 -r app && \
    useradd -u 555 -r -g app app && \
    chown -R app:app /app

# demo using non-root user

USER app

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

RUN ldd /usr/local/lib/libDRBMD.so && \
    python /app/dll/python/dll_test.py
