A1-Brouwer Degrees Working Repo
==================================

* Need to be able to do things like $(f(x)-f(y))/(x-y)$ or $f(x)/x$ when `f % ideal(x) == 0`. If we can get this to work then we can do local and global degrees over arbitrary fields
* Build trace forms along separable field extensions
* Build norms


Would be really fun:

* Implement $K_\ast^\text{M}$ and $K_\ast^{\text{MW}}$ into Macaulay2 -- by this we mean provide a nice way to give minimal presentations for elements in these rings in terms of generators, and be able to add, multiply, and check equality of elements.
* Take a form in $\text{GW}(k)$ and return the minimal $n$ so that it vanishes in $I^n$
