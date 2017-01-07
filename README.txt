This code was done for the class of G. Peyre on Sparsity and Compress Sensing.
Its provide
	- Inverse_phase: is the phase of the Fourier transform is important?
		- phase: auxiliary function
	- general_AM (Alternative Minimization) : my version of Gerchberg-Saxton
	- fourier_AM : Way faster thanks to the FFT
	- Sub_sampling : how many random intensity measurements one need to recover a signal?
	- Sparse_DR : what if we assume sparsity? (applying Douglas-Rachford)
		- DR_descent : main loop
	- Sparse_W_FB : Can we use for images with a wavelets dictionary (applying Forward-Backward)
		- FW_descent_ortho : with orthogonal wavelet
		- FW_descent : adding translation invariance to those wavelets
			*this file hasn’t been debugged*
		- vec2im : helper to handle array or vector to describe an image
	- Sparse_W_fourier_FB : what about Fourier measurements?
		- FW_descent_fourier_ortho 

Forward-Backward and Douglas-Rachford were strongly inspired by Numerical Tours available at http://www.numerical-tours.com/.
