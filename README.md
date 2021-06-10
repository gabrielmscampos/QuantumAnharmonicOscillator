# Quantum Anharmonic Oscillator

This project uses R 3.6.3 to implement a Markov Chain Monte Carlo simulation of the expected energy at the ground state of a Quantum Anharmonic Oscilator via the Virial theorem. Results were presented in the Statistical Physics class of the Rio de Janeiro State Univesity Physics Master's program.

```bash
docker run -it --rm -v $(pwd)/QuantumAnharmonicOscillator/:/mnt/qao/ r-base:3.6.3 bash

cd /mnt/qao
Rscript -e "install.packages('pracma', repos='http://cran.r-project.org')"
Rscript pathintegral.R
Rscript analyze.R
```
