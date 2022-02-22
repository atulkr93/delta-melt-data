
#!/bin/csh

autoFit.tcl -specName ft/R1rho1Dpeak-%03d.ft2 -inTab test.tab -series \
            -modX LORENTZ1D -outTab R1rho1DFit-001.tab
awk '/^set noiseRMS/ {print $4}' autoFit.com > R1rho1Dnoise-001.txt
rm autoFit.com
rm axt.tab
