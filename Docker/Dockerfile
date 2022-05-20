FROM condaforge/mambaforge
COPY environment.yml .
RUN mamba env create -f environment.yml
RUN mamba clean --all
USER 1000