# primer design notes

## prinicples of primer design paper
- codon selection - smallest change - "parsimony"
- min energy of mismatched duplex
  - mismatches treated as bulges/loops
  - total = sum(stacking energies of nearest neighbors, 
  		loops / bulges,
		dangling ends)
  - mfold 2.3 - initial nearest neighbor delta g calc
    - data with all 256 possible internal / terminal duplexes
    - in range 37-85°C - linear relationship to tm
      $ Et = a * T + b $
      T = temp
      a = slope
      b = intercept
      params a & b stored in database
- algorithm:
  1. optimal substitution
  2. design flanks - satisfy requirements in refs, similar gc content
  3. calc free energy

## general concepts paper 
- if len(primer) < 20:
	tm = 4*sum(G?C) + 2*sum(A?T)
  else:
  	tm - nearest neighbor calc


### agilent

