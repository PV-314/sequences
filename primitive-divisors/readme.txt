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
 (6) Lucas-t-2-prim-div-checks.gp
 (7) Lucas-t-3-prim-div-checks.gp
 (8) Lucas-t-4-prim-div-checks.gp
 (9) Lucas-t-6-prim-div-checks.gp

The first four files are for checking Lehmer sequences for n=5,8,10,12 respectively
(n=3,4 and 6 are not covered here, see proof in BHV)

Each of these files contains a function called
tN_check()
where N=5,8,10 or 12 is as in the name of the file.
Use these.

These tN_check() function take an optional argument, dbg.
If dbg is given a non-zero value, then extra debug information is given.

Lehmer-utils.gp contains common code used by the above four files for Lehmer sequences.

The last four files are for checking Lucas sequences for n=2,3,4,6 respectively.

Each of these files also contains a function called
tN_check()
where N=2,3,4 or 6 is as in the name of the file.
Use these.

CONTACT
=======
If you have any questions, problems, need a hand, find a bug,..., please contact PV
you can find PV's details in the arxiv preprint.
