# 1. Initialize conda in your shell (run once)
echo "source \$HOME/pl0296-02/project_data/miniforge3/etc/profile.d/conda.sh" >> ~/.bashrc
echo "module unload r openjdk 2>/dev/null" >> ~/.bashrc

# 2. Configure conda to find shared environments
cat > ~/.condarc <<'CFG'
envs_dirs:
  - $HOME/pl0296-02/project_data/aszabelska/conda-envs
  - $HOME/conda-envs
channels:
  - conda-forge
  - bioconda
auto_activate_base: false
CFG

# 3. Reload your shell config
source ~/.bashrc

# === Daily usage ===

# Activate the SUPPA2 environment
conda activate suppa2

# Verify it works
suppa.py --help
