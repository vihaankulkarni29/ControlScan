# ControlScan Epistasis Dry Run Log

Date: 2026-02-21

## Purpose
Validate the epistasis engine on a small synthetic mutation dataset to confirm co-occurrence matrix construction and Jaccard-based exclusivity scoring.

## Steps Performed
1. Created run_epistasis_test.py to exercise the MutationMatrix workflow with simulated E. coli mutation data.
2. Ran the script from the repository root to capture real output and verify the test path.

## Environment
- OS: Windows
- Python: 3.14.0 beta 2 (system install)
- Working directory: C:\Users\Vihaan\Documents\ControlScan

## Timestamps
- Script creation time: 2026-02-21 19:48:05 +05:30
- Test execution time: 2026-02-21 19:48:34 +05:30

## Command
```
C:/Users/Vihaan/AppData/Local/Programs/Python/Python314/python.exe run_epistasis_test.py 2>&1 | Out-String
```

## Output
```
Traceback (most recent call last):
  File "C:\Users\Vihaan\Documents\ControlScan\run_epistasis_test.py", line 2, in <module>
    from src.controlscan.epistasis import MutationMatrix
  File "C:\Users\Vihaan\Documents\ControlScan\src\controlscan\__init__.py", line 13, in <module>
    from .network import NetworkBuilder
  File "C:\Users\Vihaan\Documents\ControlScan\src\controlscan\network.py", line 12, in <module>
    from scipy import stats
  File "<frozen importlib._bootstrap>", line 1423, in _handle_fromlist
  File "C:\Users\Vihaan\AppData\Local\Programs\Python\Python314\Lib\site-packages\scipy\__init__.py", line 143, in __getattr__
    return _importlib.import_module(f'scipy.{name}')
  File "C:\Users\Vihaan\AppData\Local\Programs\Python\Python314\Lib\importlib\__init__.py", line 88, in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
  File "C:\Users\Vihaan\AppData\Local\Programs\Python\Python314\Lib\site-packages\scipy\stats\__init__.py", line 626, in <module>
    from ._stats_py import *
  File "C:\Users\Vihaan\AppData\Local\Programs\Python\Python314\Lib\site-packages\scipy\stats\_stats_py.py", line 52, in <module>
    from . import distributions
  File "C:\Users\Vihaan\AppData\Local\Programs\Python\Python314\Lib\site-packages\scipy\stats\distributions.py", line 16, in <module>
    from ._entropy import entropy
  File "<frozen importlib._bootstrap>", line 1371, in _find_and_load
  File "<frozen importlib._bootstrap>", line 1342, in _find_and_load_unlocked
  File "<frozen importlib._bootstrap>", line 938, in _load_unlocked
  File "<frozen importlib._bootstrap_external>", line 758, in exec_module
  File "<frozen importlib._bootstrap_external>", line 854, in get_code
  File "<frozen importlib._bootstrap_external>", line 953, in get_data
KeyboardInterrupt

Command exited with code 1
```

## Result
The dry-run test did not complete. Importing SciPy during module load was interrupted (KeyboardInterrupt). The epistasis logic did not execute because the import failed before reaching MutationMatrix runtime.

## Notes
- The test was attempted to validate epistasis behavior end-to-end from script execution.
- A follow-up run may be needed after confirming SciPy compatibility with Python 3.14 beta and allowing the import to finish.

## Stable Python 3.11 Rerun

### Purpose
Re-run the dry test using a stable Python/SciPy stack to avoid the SciPy import interruption seen with Python 3.14 beta.

### Environment
- OS: Windows
- Python: 3.11.9 (virtual environment at .venv)
- Working directory: C:\Users\Vihaan\Documents\ControlScan

### Commands
```
py -3.11 -m venv .venv
C:/Users/Vihaan/Documents/ControlScan/.venv/Scripts/python.exe -m pip install numpy==1.26.4 pandas scipy scikit-learn GEOparse
C:/Users/Vihaan/Documents/ControlScan/.venv/Scripts/python.exe run_epistasis_test.py 2>&1 | Out-String
```

### Timestamps
- Environment setup and dependency install: 2026-02-21 20:02:36 +05:30
- Final test run: 2026-02-21 20:03:40 +05:30

### Interim Failures Observed (and Resolved)
- 2026-02-21 19:54:39 +05:30: NumPy C-extension mismatch (cp314 artifacts in a cp311 environment). Fixed by recreating the venv and reinstalling with cp311 wheels.
- 2026-02-21 19:59:54 +05:30: pandas missing because the first install attempt was cancelled. Re-ran install to completion.
- 2026-02-21 20:03:40 +05:30: Updated MutationMatrix to accept lowercase input keys ('mutations') to match the test data format.

### Output
```
--- ControlScan Epistasis Engine Test ---

[1] Building Co-Occurrence Matrix...
            acrA_T104A  acrB_N596H  gyrA_S83L  parC_S80I
acrA_T104A           5           4          1          1
acrB_N596H           4           4          1          1
gyrA_S83L            1           1          2          1
parC_S80I            1           1          1          2


[2] Calculating Epistatic Probabilities...
Jaccard Co-occurrence Score for acrA_T104A & acrB_N596H: 0.80
>> BIOLOGICAL CONCLUSION: High Epistasis Detected. These mutations are structurally co-dependent.
```

### Result
The stable Python 3.11 rerun completed successfully and produced the expected co-occurrence matrix and epistasis score.
