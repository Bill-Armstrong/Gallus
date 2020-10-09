# Gallus
## Efficient, understandable and fast machine learning.

This is a challenge project for persons familiar with C++ who want to learn how to program a GPU using CUDA to do machine learning. The task uses the repository NANO (Noise-Attenuating Neuron Online):

git clone https://github.com/Bill-Armstrong/NANO.git

The repository contains an ALN machine learning program. ALN technology is different from the current neural network technology. An ALN consists of a first layer computing scalar products of an input vector with weight vectors of the nodes. All other layers are maximum and minimum operators. The learned functions are piecewise linear and continuous.  It is able to do classification of patterns captured on a retina. For example it has been used (without CUDA acceleration) to classify handwritten digits in the well-known MNIST data set to very high accuracy. For a given target class it creates a membership function that slowly decreases as the L1 distance from samples of the target class increases. With some minor changes the NANO project is expected to be able to store a large, growing number of image classes, and be able to classify new images quickly in almost constant time and with a constant amount of hardware except for storage. The goal of GALLUS is to accelerate the NANO project code to speed up the scalar products in the first layer of the neural network. HINT: You might adapt an NVIDIA CUDA Toolkit example called scalarProd which doeas most of the work.

The Gallus project originated in a university laboratory to train students in machine learning.  It was then called "Pullus: the electric chicken". "Pullus", in Latin, refers to a young chick in the downy stage. The project is now full-grown and "Gallus" is the Latin word for rooster. The 1953 drawing by Bernard Buffet shows the rooster as proud and indomitable!
