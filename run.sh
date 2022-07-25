nextflow kuberun xsvato01/mareckova_CXCR_k8s -pod-image 'cerit.io/nextflow/nextflow:22.05.0'\
	-c zaloha_nextflow.config -v pvc-janek-storage-elixir1-cerit-sc-cz:/mnt \
        --datain /mnt/home/450402/000000-My_Documents/mareckova_data/ --outdir /mnt/home/450402/000000-My_Documents/mareckova_results/
