FROM condaforge/mambaforge
COPY environment.yml .
RUN mamba env create -f environment.yml
RUN conda activate mareckova_CXCR_k8s
USER 1000