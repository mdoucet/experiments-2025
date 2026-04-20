# Copper film

## Description
Deposited 50 nm copper on 3 nm titanium on silicon, in D2O.
Incoming beam is coming from the back of sample (back reflection).
Copper oxide is likely present on the copper surface.

## Details
- Electrolyte: D2O + 0.1 M NaHCO3 sparged with N2
- Potential experimental issue: Electrolyte was sparged for too long (~10 min?), couldn't find pH probe or pH paper so took out a small aliquot and put in a sealed vial. 
- pH measured via pH paper = 9 
- Alkaline conditions could cause persistent hydroxyls on the surface.
- Sample was aligned once at the beginning of the experiment.

## Measurements:
- Run 226642: OCV
- Run 226645: -0.9 V vs. Ag/AgCl tNR 
- Run 226646: -0.9 V vs. Ag/AgCl tNR 
- Run 226649: -0.9 V vs. Ag/AgCl tNR 
- Run 226652: Return to OCV tNR
- Run 226553: OCV
- Run 226656: OCV
- Run 226659: OCV 

# Fitting approach
- Co-refine segments, not the combined data sets.
- Allow for sample broadening.
- Allow for angle offset (same for all segments of a measurement).

# Fits to perform 
1. Run 226642
    - Initial OCV
2. Run 226649
    - Time resolved data that needs to be sliced in 30-sec time intervals.
3. Corefine runs 226642 and 226649
    - Keep parameters tied except the parameters related to the copper oxide and the titanium.
    - The solvent SLD, and the copper SLD should be the same for both measurements.
    - The thickness of the copper may vary slightly.

## Questions
- Are outcomes of fits 1 and 2 in agreement with fit 3?
- If I'm right about pH = 9 inducing oxidation..
    1. does oxide reduce during tNR at all?? 
    2. does oxide grow or does surface roughen in 226553-226659? 