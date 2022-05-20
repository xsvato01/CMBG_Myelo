nextflow kuberun xsvato01/nextflow_k8s_test/my_pipeline.nf -pod-image 'cerit.io/nextflow:21.09.1' \
        -w /mnt/home/450402/000000-My_Documents/k8s_nf/tmp \
        -v pvc-janek-storage-elixir1-cerit-sc-cz:/mnt -c zaloha_nextflow.config \
        --outdir /mnt/home/450402/000000-My_Documents/Next_ik refindex  /mnt/shared/999993-Bioda/data/ssd_3/references/homsap/GRCh38-p10/index/BWA/GRCh38-p10 \
        --refvcf /mnt/shared/999993-Bioda/data/ssd_3/references/homsap/GRCh38-p10/seq/GRCh38-p10.fa