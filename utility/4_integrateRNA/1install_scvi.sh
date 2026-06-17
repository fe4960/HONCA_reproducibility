source ~/.condainit
conda clean --all --force-pkgs-dirs --yes
mamba create -n scvi python=3.10 bash pytorch-gpu parallel leidenalg blas=*=openblas anndata scanpy matplotlib pandas scvi-tools 'pulp<2.8'
conda activate scvi
download https://github.com/RCHENLAB/SingleCellAnalysis/tree/main
pip install /pub/junw42/software/SingleCellAnalysis-main/Modules/Bash/jlbashwrapperpubutils-0.0.1-py3-none-any.whl
source $(jlbashwrapperpubutils -e)
pip install /pub/junw42/software/SingleCellAnalysis-main/Modules/Python/jlpythonwrapperpubutils-0.0.1-py3-none-any.whl
