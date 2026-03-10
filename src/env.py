import os
from pathlib import Path

def _get_env(name: str) -> str:
    value = os.getenv(name)
    if value is None:
        raise EnvironmentError(f'Environment variable {name} not set')
    return value

DATA_ROOT = Path(_get_env('DATA_ROOT'))
