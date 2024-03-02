# HBStools C algorithms

In this folder there you'll find source code to a library implementing:

* A pure implementation of Poisson-FOCuS, with all optimizations we are aware of. To do its thing, the algorithm expect you to provide observed counts and expected background rates.
* Poisson-FOCuS+SES, a wrapper to Poisson-FOCuS, providing automatic background estimate via simple exponential smoothing. This will require you to only provide observed counts.
* The BFT, a wrapper to an arbitrary number of Poisson-FOCuS+SES algorithms. The role of the BFT is to orchestrate these algorithms together. It requires you providing a number of counts, one per each sub-algorithm.

They are put together as a matrioska. The BFT contains many PF+SES, and each PF+SES contains a PF algorithm.
BFT, PF+SES, and PF all have the same interfaces, these "matrioska" are equally painted.
The approach to runtime errors is that there should be none. If an algorithm encounter one, it tells you, and do not perform any dangerous operations. If you keep providing shit data, the algorithm keeps telling you and do nothing.
This is because these implementations were initially designed for safety critical applications. 

![bft](../../../assets/matrioska.jpeg)


A few more words on this. The actual interface (the "offline" API) to this library used by HBStools returns errors if enough shit data are passed through. 
This is because HBStools is intended to be used offline, so there are no dangers involved with errors.
The underlying code (the "online API", which is the interface intended for safety-critical applications) implements the no-runtime-errors approach talked above.
So can this code, and its online API, be copy-pasted on whatever embedded microcontroller orbiting Earth or on your car without you dying? 
In practice yes, but keep in mind we are using dynamic allocation for memory heap chunks of constant size. Again, we do this because the most dangerous thing this implementation will encounter is the coffee cup on your desktop: again, HBStools runs offline. However, the code is designed to make transition to a compile-time allocation trivial, defining your own memory pool and replacing the mallocs with some addressing to the pool. 
