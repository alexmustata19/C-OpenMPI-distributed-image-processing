build: filtru

filtru: tema3.c
	mpicc tema3.c -lm -o filtru

clean:
	rm -f filtru