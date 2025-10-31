conda activate nfc

cd /Users/daisy/home/sviluppo/hackathon-oct-2025

git clone git@github.com:nf-core/sammyseq.git

cd /Users/daisy/home/sviluppo/hackathon-oct-2025/sammyseq

git checkout rustem_fixes
git checkout add_compartments
git checkout subworkflow_update


cd /Users/daisy/home/sviluppo/hackathon-oct-2025
nextflow run sammyseq -profile test --outdir test_before_commit
rsync -av sammyseq mmutarelli@nvidia.hpc.isasi.:sviluppo/

cd /Users/daisy/home/sviluppo/hackathon-oct-2025/sammyseq
pre-commit run --all-files
nf-core pipelines lint
nf-core pipelines lint > ../lint_report.txt
code ../lint_report.txt



ssh mmutarelli@nvidia.hpc.isasi.

cd ~/sviluppo
nextflow run sammyseq -profile test --outdir test_before_commit --comparison S2vsS3
