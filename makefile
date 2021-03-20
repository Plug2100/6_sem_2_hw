test: main.cpp
	mpic++ -o main main.cpp -lm
	for n in 25 26 27; do \
    		for k in 1 13 $$n; do \
    			for i in 1 2 4 8 16 32 64; do \
    				mpirun -np $$i main matr.bin k $$k n $$n file_write outp.bin ; \
    			done \
    		done \
    	done
    	
    	
    	
    	
check: main.cpp
	mpic++ -o main main.cpp -lm
	mpirun -np 1 main k 1 n 2  file_read inp.bin test file_test check.bin
	mpirun -np 2 main k 1 n 2  file_read inp.bin test file_test check.bin
	mpirun -np 4 main k 1 n 2  file_read inp.bin test file_test check.bin
	mpirun -np 8 main k 1 n 2  file_read inp.bin test file_test check.bin
	mpirun -np 16 main k 1 n 2  file_read inp.bin test file_test check.bin
