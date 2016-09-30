# GreenBST: Energy-efficient concurrent search tree

This repository contains the source code and the benchmark framework from the below article:

`
Ibrahim Umar, Otto J. Anshus, and Phuong H. Ha. GreenBST: An energy-efficient concurrent search tree. Proceedings of the 22nd International European Conference on Parallel and Distributed Computing (Euro-Par ’16), 2016, LNCS, pp. 502-517, Springer.[(link)](http://dx.doi.org/10.1007/978-3-319-43659-3_37)
`

### License

* GreenBST, the benchmark routines, and the building frameworks is developed by UiT and licensed under the Apache License, Version 2.0.

* CBTree and SVEB are developed by UiT based on the others’ work. They are licensed under the The GNU General Public License v3.0.

* Other trees are developed by their respective authors and retain their original licenses

### To run the benchmark:

1. Please make sure you have the Intel PCM library and PAPI library installed.
2. Edit the ```common/common.mk``` to modify the the Intel PCM and PAPI libraries location.
3. Go to ```bench/``` directory.
4. Run ```make-bins.sh```. This should create the trees' binaries.
5. Run ```runtest.{ARCH}.sh``` to run the benchmark.
6. You will find the benchmark results (in CSV format) inside the ```bench/results/``` directory

### The trees

**1. GreenBST**

GreenBST is a locality-aware and energy-efficient concurrent search tree. GreenBST portable, namely it can maintain their locality-awareness and energy-efficiency on different computing platforms (platform-independent).

**2. Concurrent B-tree (CBTree)**

CBTree is a prominent locality-aware concurrent B+tree. CBTree is a representation of the classic coarse-grained locality-aware search concurrent trees that are usually platform-dependent. CBTree can only perform well if their node size is set correctly (e.g., to the system’s page size). This tree is also often referred as the B-link tree.

* Philip L. Lehman and s. Bing Yao. 1981. Efficient locking for concurrent operations on B-trees. ACM Trans. Database Syst. 6, 4 (December 1981), 650-670.

**3. Lock-based (SVEB) dynamic cache-oblivious tree**

SVEB is the concurrent implementation of the fine-grained locality-aware vEB binary search tree. SVEB uses a global mutex to serialize its concurrent tree operations.

* Gerth Stølting Brodal, Rolf Fagerberg, and Riko Jacob. Cache oblivious search trees via binary trees of small height. In Proc. 13th ACM-SIAM Symp. Discrete algorithms, SODA ’02, pages 39–48, 2002.

**4. Non-blocking binary search tree (LFBST)**

LFBST is the improved variant of the original non-blocking binary search tree. These non-blocking search trees are locality-oblivious.

*LFBST official repository:*
https://github.com/anataraja/lfbst

* Aravind Natarajan and Neeraj Mittal. Fast concurrent lock-free binary search trees. In Proc. 19th ACM SIGPLAN Symposium on Principles and Practice of Parallel Programming, PPoPP ’14, pages 317–328, 2014.

**6. RCU-based concurrent search tree (Citrus)**

Citrus is a concurrent binary search tree that utilizes Read-Copy-Update (RCU) synchronization and fine-grained locking for synchronization among updaters. Citrus contain operation is wait-free. This concurrent search tree is locality-oblivious.

*Citrus official repository:*
https://bitbucket.org/mayaarl/citrus

* Maya Arbel and Hagit Attiya. Concurrent updates with rcu: Search tree as an example. In Proc. 2014 ACM Symposium on Principles of Distributed Computing, PODC ’14, pages 196–205. ACM, 2014.

**7. Portably scalable concurrent search tree (BSTTK)**

BST-TK is the state-of- the-art lock-based concurrent search tree based on the asynchronous concurrency paradigm. BST-TK is portably scalable, namely it scales across different types of hardware platforms. BSTTK is a locality-oblivious tree.

*BSTTK official repository:*
https://github.com/LPD-EPFL/ASCYLIB

* Tudor David, Rachid Guerraoui, and Vasileios Trigonakis. Asynchronized concurrency: The secret to scaling concurrent search data structures. In Proc. 12th Intl. Conf. on Architectural Support for Programming Languages and Operating Systems, ASPLOS’15, pages 631–644, 2015
