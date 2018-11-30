# Pull base image.
FROM brmather/pybadlands-dependencies:0.1

MAINTAINER Ben Mather

RUN apt-get remove -y python3-six
RUN python3 -m pip install --upgrade six requests html5lib urllib3 matplotlib

WORKDIR /build
RUN git clone https://github.com/brmather/pyBadlands.git && \
    cd pyBadlands && \
    git checkout python3
WORKDIR /build/pyBadlands/pyBadlands/libUtils
RUN make
RUN python3 -m pip install -e /build/pyBadlands

RUN git clone https://github.com/matplotlib/legacycontour.git && \
    python3 -m pip install -e /build/legacycontour

WORKDIR /build
RUN git clone https://github.com/badlands-model/pyBadlands-Companion.git
RUN python3 -m pip install -e /build/pyBadlands-Companion

ENV TINI_VERSION v0.8.4
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/local/bin/tini
RUN chmod +x /usr/local/bin/tini

# Copy cluster configuration
RUN mkdir /root/.ipython
COPY profile_mpi /root/.ipython/profile_mpi

RUN mkdir /workspace && \
    mkdir /workspace/volume && \
    mkdir /workspace/companion

# Copy test files to workspace
RUN cp -av /build/pyBadlands/Examples/* /workspace/
RUN cp -av /build/pyBadlands-Companion/notebooks/* /workspace/companion/

COPY run.sh /build
RUN chmod +x /build/run.sh

# setup space for working in
VOLUME /workspace/volume
# launch notebook
WORKDIR /workspace
EXPOSE 8888
ENTRYPOINT ["/usr/local/bin/tini", "--"]

RUN find -name \*.ipynb -print0 | xargs -0 jupyter trust

ENV LD_LIBRARY_PATH=/workspace/volume/pyBadlands/pyBadlands/libUtils:/build/pyBadlands/pyBadlands/libUtils
CMD /build/run.sh