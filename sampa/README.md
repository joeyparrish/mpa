# Original SAMPA

This is (nearly) the original SAMPA (Solar and Moon Position Algorithm) code as
downloaded from the US Department of Energy and developed by the National
Renewable Energy Laboratory (NREL).

We use it here for regression testing, so that we are certain our changes do
not impact accuracy of results.

In order to facilitate comparison testing, this code has been modified in the
following ways:

 - .c files renamed to .cc to allow namespacing
 - "namespace sampa" added to avoid collisions with the new implementation when
   both are linked into the regression test executable
 - Automated formatting applied for readability

**TODO**: Add a link to the original download location from the DoE
