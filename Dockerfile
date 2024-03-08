FROM continuumio/miniconda3:23.5.2-0

RUN apt-get -y update && apt-get -y install gcc libcfitsio-bin && \
        apt-get autoclean && \
        rm -rf /var/lib/apt/lists/*
        
RUN conda install -y numpy astropy pytest pytest-runner ipython matplotlib scipy \
        && conda install watchdog emcee corner -c conda-forge \
        && conda clean -y --all

RUN mkdir /home/eng \
        && /usr/sbin/groupadd -g 500 "eng" \
        && /usr/sbin/useradd -g 500 -d /home/eng -M -N -u 500 eng \
        && chown -R eng:eng /home/eng

RUN pip install lcogt-logging Cython \
        && rm -rf ~/.cache/pip

RUN git clone https://github.com/mstamy2/PyPDF2 /usr/src/pypdf2
RUN pip install /usr/src/pypdf2

RUN git clone https://github.com/kbarbary/sep.git /usr/src/sep
RUN pip install /usr/src/sep

RUN git clone https://github.com/astropy/astroscrappy.git /usr/src/astroscrappy
RUN pip install /usr/src/astroscrappy

COPY . /pupepat/src/

RUN pip install /pupepat/src

ENV HOME /home/eng

WORKDIR /home/eng

USER eng
