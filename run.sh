nextflow kuberun xsvato01/mareckova_CXCR_k8s/ -pod-image 'cerit.io/nextflow:21.09.1' \
        -w /mnt/home/450402/000000-My_Documents/mareckova_k8s_work/tmp \
        -v pvc-janek-storage-elixir1-cerit-sc-cz:/mnt -c nextflow.config \
        --datain /mnt/home/450402/000000-My_Documents/mareckova_data/