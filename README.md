# Baluyot-Calibration-Algorithm-Implementation
This repository contains a simple implementation of the Baluyot Calibration Algorithm on C++.

This implementation was used for obtaining experimental results to prove the algorithm's correctness and accompanies the article on:
(https://medium.com/@bastien.baluyot/the-baluyot-calibration-algorithm-85c16bc1c878)

This implementation is also based on the original context of measuring a battery's SoC (though this term is only to represent the input and is actually unbounded) and assumes 2 inputs (SoC1 and SoC2), i.e. n = 2, as explained in the article.

The code allows for 3 input settings: Manual (via cin input), Automated (via PRNG) or Automated (via external .txt file parse).
All relevant output data will be logged and stored in external .txt files.
