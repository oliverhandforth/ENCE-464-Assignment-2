all: threads

# -g outputs debugging information
# -pthread adds threading support
threads: threads.c
	gcc -g -pthread -o threads threads.c

clean:
	rm -f threads *.o
