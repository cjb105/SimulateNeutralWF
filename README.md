# SimulateNeutralWF
Programs to accompany the paper "An approximate solution for multi-allele neutral diffusion with low mutation rates"

R programs to simulate the neutral multi-allele Wright-Fisher model:

simulateNeutralWF3Alleles.R produces Figure 3;

simulateNeutralWF4Alleles.R produces Figures 5 and 6;

The parameters 
  N (population size), 
  alpha12 etc., Pi and phi (parameterising the rate matrix)
can be adjusted, but not K (number of alleles)

Mathematica source codes to derive Eqs. (27) and (46), that is, the result
  C_ab = 2 alpha_ab pi_a pi_b 
for the normalisation of effective line densities on edges of the solution simplex:

Ccompute-3alleles.nb           does the K=3 alleles case;

Ccompute-4alleles-leaveout1.nb does the K=4 alleles case
