Since we are dealing with paired end data, whenever both reads from the pair are mapped, they will be reported as separate entries in the BAM file, but they will share the read name. 
On the other hand, if only one read of the pair is mapped, it will be reported only once.
