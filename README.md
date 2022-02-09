# casanger

Cas9 nickase validation by Sanger sequencing analysis workflow: CaSanger.

## Background

CRISPR Cas9 is a RNA guided endonuclease now commonly utilized in a broad range
of molecular biology applications. The native Cas9-sgRNA complex produces DNA
double strand breaks when successfully targeted. Variants of Cas9 have
been engineered with only one active catalytic site; resulting in th enzyme
producing single strand DNA breaks (nicks) when successfully targeted. These
Cas9 "nickases" are commonly used with 2 sets of sgRNAs which together
produce double strand breaks. 

A key component of any Cas9 assay is validation of both the presence and location
of a specific break. There are a variety of methods to validate the presence
of double strand DNA breaks both *in vivo* and *in vitro* but validation of
a nick is more difficult largely due to difficulty of designing a PCR reaction
that will produce target, and nick dependent products.

An alternative approach is to utilize Sanger sequencing which is sensitive 
to DNA perturbations and relies on linear amplification (utilizing a single
primer). If a nicked DNA template is subjected to Sanger sequencing with a
primer that binds to the same strand that is targeted by a nicking
enzyme there should be no distinguishable signal beyond the nick. This
is because if a sufficiently high percentage of the template has the same
targeted nick, the polymerase used in the sequencing reaction will fall off
the single stranded template at the site of nick.

## Running the workflow


### Configuration files


## Output






