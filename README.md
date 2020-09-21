# Gallus
More efficient and understandable machine learning with C++ and CUDA

IMPORTANT: THIS PROJECT DOES NOT WORK YET BECAUSE OF A LINKAGE PROBLEM  BETWEEN C++ and CUDA.

When it does work, it will be able to do classification of sensed patterns. For example it has been used (without CUDA acceleration) to classify handwritten digits in the well-known MNIST data set to very high accuracy. For a given target class it creates a membership function that slowly decreases as the L1 distance from samples of the target class increases. When it is finally working and complete, it is expected to be able to store a large, growing number of image classes, and be able to classify new images quickly in almost constant time and with a constant amount of hardware except for storage.

The project originated in a university laboratory to train students in machine learning.  It was then called "Pullus: the electric chicken". "Pullus", in Latin, refers to a young chick in the downy stage. The project is now full-grown and "Gallus" is the Latin word for rooster. The 1953 drawing by Bernard Buffet shows the rooster as proud and indomitable!
