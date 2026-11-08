"""One-command reproduction of the benchmark and all quality controls."""
from pathlib import Path
import subprocess
import sys

P=Path(__file__).resolve().parent
for script in ('solution.py','verification.py','test_physics.py','make_plots.py'):
    print(f'\nRUN {script}',flush=True)
    subprocess.run([sys.executable,str(P/script)],cwd=P,check=True)
print('\nCOMPLETE: All configured runs and tests succeeded.',flush=True)
