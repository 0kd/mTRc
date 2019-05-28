# Makefile

PROGRAM = mTRc
OBJS	= main.o k_means_clustering.o print_and_feed_one_repeat.o 
CC	= cc
CFLAGS	=

.c.o:
	$(CC) $(CFLAGS) -c $<

$(PROGRAM): $(OBJS)
	$(CC) $(OBJS) -o $(PROGRAM)

clean:
	rm $(PROGRAM) $(OBJS)