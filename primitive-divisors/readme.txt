INSTALLATION
============
(1) have Pari installed
(2) in the directory where Pari is installed, create a directory called sequences with a subdirectory called primitive-divisors
(3) put all the files in https://github.com/PV-314/sequences/primitive-divisors in the primitive-divisors directory

CONTENTS
========
 (1) Lehmer-t-5-prim-div-checks.gp
 (2) Lehmer-t-8-prim-div-checks.gp
 (3) Lehmer-t-10-prim-div-checks.gp
 (4) Lehmer-t-12-prim-div-checks.gp
 (5) Lehmer-utils.gp

The first four files are for checking n=5,8,10,12 respectively
(n=3,4 and 6 are not covered here, see proof in BHV)

Each of these files contains a function called
tN_check()
where N=5,8,10 or 12 is as in the name of the file.

These tN_check() function take an optional argument, dbg.
If dbg is given a non-zero value, then extra debug information is given.

Finally, Lehmer-utils.gp contains common code used by the other four files.

CONTACT
=======
If you have any questions, problems, need a hand, find a bug,..., please contact PV
you can find PV's details in the arxiv preprint.
