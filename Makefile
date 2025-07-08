CC			:= gcc
CFLAGS		:= -O3 -lm -I. -DCQT_DEBUG -DSM_DEBUG

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: pffft/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: test_cqt.o pffft.o cqt.o
	$(CC) $(CFLAGS) $^ -o $@
	./$@

visualize: test_visualize.o pffft.o cqt.o
	$(CC) $(CFLAGS) $^ -o $@
	./$@

.PHONY: clean
clean:
	rm -r *.o test visualize
