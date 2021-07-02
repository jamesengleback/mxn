# `mxn`

Primer design for site-directed mutagenesis

## background

Site-directed mutagenesis is a relatively cheap way of making point mutations in a DNA sequence (if it works). It works by PCR of a circular plasmid containing the gene of interest where the mutations are carried by the primers and should be present in all PCR copies of the gene. 

Manual primer design is time consuming and doesn't scale well. Online tools are ok but don't integrate well with other tools. `mxn` is an attempt to automate primer design for this purpose, which can allow for more complex reactions that aren't limited by human design.

### PCR methods

Different site-directed mutagenesis kit providers have different primer design  requirements and those primers may have a different melting temperature (Tm) depending on the contents of the reaction kit (e.g. DMSO lowers the effective Tm). The two kits that `mxn` was designed for are *NEB Q5* and *Agilent QuickChange Multi*.

#### *NEB Q5*
*NEB Q5* kits use back to back primers where the mutation payload is carried in the tail of one. The full sequence is copied from circular to linear double strands with the new mutation at one end. The linear DNA is recircularized using a Kinase and a Ligase and the original copy is digested by a DPN1. 

TM requirements, calculator & attempts to replicated, primers

#### *Agilent QuickChange Multi*

*Agilent QuickChange Multi* kits can introduce up to about 4 mutations in one reaction which is cool. This is a single strand PCR where the mutations are carried mid-primer. Nicks are sealed apparently and the product is sing strand circular DNA, which you can transform.

TM requirements, calculator & attempts to replicated, primers


## Tm estimation

## features

## howto

