CC = gcc
CFLAGS = -mcmodel=medium -Wall -O2 -std=c99 -fpic 

OBJ = SynSearch.o 

synsearch: $(OBJ)
	$(CC) $(CFLAGS) -o synsearch SynSearch.o

clean:
	rm -f *.o synsearch


