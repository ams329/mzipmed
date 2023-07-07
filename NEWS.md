#mzipmed 1.4.0
-Added offset options to mzip and mediation functions
-Function with binary outcome and ZI mediator can now use logistic regression outcome model using the `rare` argument
-Adding an offset added made some derivations conditional on a fixed value of the offset. These functions had an additional argument called `OFF` added
-Removed dependence on `robust` package and added dependence on `sandwich` package to be used when computing the robust covariance matrix for the Poisson outcome model

---

#mzipmed 1.3.5
-Fixed typos in the vignette 
-Added some additional information to the vignette for clarity purposes
-Other minor code fixes

---

#mzipmed 1.3.3 & 1.3.4
-Various typo fixes

---

#mzipmed 1.3.2
-Corrected some terminology issues relating to risk ratios vs. rate ratios

---

#mzipmed 1.3.1
-Made 'print' default in `mzip` function TRUE
-Fixed print in `mzip` function to not give 'relative risk' of intercept

---

# mzipmed 1.3.0
-Added vignette

---
# mzipmed 1.2.7
-Added R console output for `zioutlmmedint` that was accidentally left out prior

---
# mzipmed 1.2.6
-Added R console output for mediation functions and beautified for MZIP function

---

# mzipmed 1.2.5

-Made some edits to the description file to satisfy CRAN requirements
-Various grammatical fixes to documentation

---

# mzipmed 1.2.4

-Added contact information of package maintainer

---
# mzipmed 1.2.3

-Removed package dependency 'boot'
-Added citation for MZIP manuscript

---

# mzipmed 1.2.2

-Fixed some issues identified through R CMD check
-Changed Leann Long from contributor to author

---

# mzipmed 1.2.1

-Added a `NEWS.md` file to track changes to the package.
-Added MIT license
-Fixed some naming/character/spacing issues
-Delta method mediation can now use MZIP with robust covariance matrix
-MZIP function now outputs robust covariance matrix
-Added a simulated data file called mzipmed_data
-Changed example to be useable with mzipmed_data


---

# mzipmed 1.2.0

-Added functions for mediation with zero-inflated outcomes using MZIP
-Fix typos in documentation


---

# mzipmed 1.1.1

-Deleted some unnecessary lines from mzip function to make computation quicker

---

# mzipmed 1.1.0

-Fixed some formulas for binary/count outcome after discovering error in derivations