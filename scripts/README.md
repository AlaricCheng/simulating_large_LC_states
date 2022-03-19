- `lib/` is the directory for necessary libraries. 


- `fusion_circuits.m` converts 1200 raw records to probability distributions which are used in the following python scripts.  
Add `lib` to matlab path then run the script in the `scripts` directory.


- `likelyrho_cluster12.m` processes the probability distributions from 12-qubit linear cluster states with the maximum likelihood method.  
Add `lib` to matlab path then run the script in the `scripts` directory. Note that it will take about 8 min.


- `fig_subcircuit_data.py` draws Fig. 3 and Fig. S1 (a). Run 
```bash
python ./scripts/fig_subcircuit_data.py
```
in the root directory. 


- `fig_large_LC_data.py` draws Fig. 4 and Fig. S1 (b). Run 
```bash
python ./scripts/fig_large_LC_data.py
```
in the root directory. Note that it will take about 5 min.


- `fig_waveform_1.py` draws Fig. 2 (b). Run 
```bash
python ./fig_waveform_1.py
```
in the `scripts` directory.


- `fig_waveform_2.py` draws Fig. 2 (c). Run 
```bash
python ./fig_waveform_2.py
```
in the `scripts` directory.


- `stab_calibration.py` provides the text-form-data. Run 
```bash
python ./stab_calibration.py
```
in the `scripts` directory.