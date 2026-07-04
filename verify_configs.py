from pathlib import Path
from dsvr.config import load_config

config_dir = Path("configs")
for yaml_file in config_dir.glob("*.yaml"):
    if yaml_file.name == "auto3d_entropy.auto3d.yaml":
        continue
    print(f"Loading {yaml_file}...")
    try:
        config = load_config(yaml_file)
        print(f"  Description: {config.description}")
    except Exception as e:
        print(f"  FAILED to load {yaml_file}: {e}")
        exit(1)
print("All configs loaded successfully!")
